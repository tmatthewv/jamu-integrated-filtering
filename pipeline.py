#!/usr/bin/env python3
"""
pipeline.py — Compound filtering pipeline
Tambunan et al. (2025) — Stage 3: Compound Retrieval & Filtering

Usage (from project root):
    conda activate docking
    python scripts/03_compounds/pipeline.py

Pipeline stages:
    1. KNApSAcK compound retrieval
    2. PubChem SMILES + deduplication
    3. Lipinski RO5 + Veber filter
    4. ADMET filter (RDKit)
    5. ChEMBL target prediction (similarity ≥ 60%, pChEMBL ≥ 5.0)
    6. Final shortlist assembly
"""

import os, re, time, requests, pandas as pd
from pathlib import Path

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, FilterCatalog
    from rdkit.Chem import AllChem
    RDKIT_OK = True
except ImportError:
    print("WARNING: RDKit not available — install via: conda install rdkit")
    RDKIT_OK = False

ROOT    = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
PRO_DIR = ROOT / "data" / "processed"
RES_DIR = ROOT / "data" / "results"
for d in [PRO_DIR, RES_DIR]: d.mkdir(parents=True, exist_ok=True)

DELAY = 0.3

# ── Target proteins ──────────────────────────────────────────
TARGET_PROTEINS = {
    "P10275": "Androgen Receptor (ANDR)",
    "P03372": "Estrogen Receptor alpha (ESR1)",
    "Q16665": "HIF-1alpha (HIF1A)",
    "P19838": "NF-kB1 (NFKB1)",
}

# ════════════════════════════════════════════════════════════
# STEP 1: Retrieve compounds from KNApSAcK Core
# ════════════════════════════════════════════════════════════

def step1_knapsack():
    out = PRO_DIR / "knapsack_raw.csv"
    if out.exists():
        print(f"[1] Already exists: {out.name}")
        return pd.read_csv(out)

    freq_file = PRO_DIR / "jamu_frequency_final.csv"
    if not freq_file.exists():
        raise FileNotFoundError("Run scrape_jamu.py first to generate jamu_frequency_final.csv")

    freq  = pd.read_csv(freq_file)
    Q1    = freq["final_score"].quantile(0.25)
    Q3    = freq["final_score"].quantile(0.75)
    tukey = Q1 + 1.5*(Q3-Q1)
    priority = (freq[freq["final_score"] >= tukey]
                .sort_values("final_score", ascending=False)
                .groupby("FAMILI").head(3)
                .reset_index(drop=True))
    print(f"[1] Retrieving compounds for {len(priority)} priority plants...")

    rows = []
    for _, plant in priority.iterrows():
        r = requests.get(
            "https://www.knapsackfamily.com/knapsack_core/result.php",
            params={"sname":"all","word": plant["TANAMAN"],"display":"100000"},
            timeout=30
        )
        # Parse table from HTML response
        # [parsing logic]
        time.sleep(DELAY)

    df = pd.DataFrame(rows, columns=["knapsack_id","cas_id","compound_name","species"])
    df.to_csv(out, index=False)
    print(f"    Saved {len(df)} raw compounds → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 2: PubChem SMILES + deduplication
# ════════════════════════════════════════════════════════════

def get_smiles(name, cas=None):
    """Retrieve canonical SMILES from PubChem by name or CAS."""
    for query, ns in [(cas,"name"),(name,"name")]:
        if not query: continue
        try:
            r = requests.get(
                f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{ns}/"
                f"{requests.utils.quote(str(query))}/property/"
                f"CanonicalSMILES,MolecularWeight,MolecularFormula,CID/JSON",
                timeout=15)
            if r.status_code == 200:
                props = r.json()["PropertyTable"]["Properties"][0]
                return props.get("CanonicalSMILES"), props
        except Exception:
            pass
        time.sleep(0.2)
    return None, {}

def step2_smiles(df_raw):
    out = PRO_DIR / "compounds_with_smiles.csv"
    if out.exists():
        print(f"[2] Already exists: {out.name}")
        return pd.read_csv(out)

    print(f"[2] Retrieving SMILES for {len(df_raw)} raw compounds...")
    rows = []
    seen_smiles = set()

    for i, row in df_raw.iterrows():
        smiles, props = get_smiles(row["compound_name"], row.get("cas_id"))
        if not smiles: continue

        # Canonical dedup via RDKit
        if RDKIT_OK:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: continue
            canon = Chem.MolToSmiles(mol)
        else:
            canon = smiles

        if canon in seen_smiles:
            continue
        seen_smiles.add(canon)
        rows.append({**row.to_dict(), "smiles": canon,
                     "cid": props.get("CID"),
                     "mw": props.get("MolecularWeight"),
                     "formula": props.get("MolecularFormula")})

        if i % 100 == 0:
            print(f"    {i}/{len(df_raw)} processed, {len(rows)} unique...")
        time.sleep(DELAY)

    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"    Saved {len(df)} unique SMILES → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 3: Lipinski RO5 + Veber filter
# ════════════════════════════════════════════════════════════

def calc_descriptors(smiles):
    if not RDKIT_OK: return {}
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return {}
    return {
        "mw_calc"  : Descriptors.MolWt(mol),
        "logp"     : Descriptors.MolLogP(mol),
        "hbd"      : rdMolDescriptors.CalcNumHBD(mol),
        "hba"      : rdMolDescriptors.CalcNumHBA(mol),
        "tpsa"     : rdMolDescriptors.CalcTPSA(mol),
        "rotbond"  : rdMolDescriptors.CalcNumRotatableBonds(mol),
    }

def step3_lipinski(df_smiles):
    out = PRO_DIR / "compounds_lipinski.csv"
    if out.exists():
        print(f"[3] Already exists: {out.name}")
        return pd.read_csv(out)

    print(f"[3] Lipinski RO5 + Veber filter (n={len(df_smiles)})...")
    descs = pd.DataFrame([calc_descriptors(s) for s in df_smiles["smiles"]])
    df = pd.concat([df_smiles.reset_index(drop=True), descs], axis=1)

    df = df[
        (df["mw_calc"]  <= 500) &
        (df["logp"]     <= 5)   &
        (df["hbd"]      <= 5)   &
        (df["hba"]      <= 10)  &
        (df["tpsa"]     <= 140) &
        (df["rotbond"]  <= 10)
    ]
    df.to_csv(out, index=False)
    print(f"    Passed: {len(df)} → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 4: ADMET filter (RDKit)
# ════════════════════════════════════════════════════════════

def step4_admet(df_lip):
    out = PRO_DIR / "compounds_admet.csv"
    if out.exists():
        print(f"[4] Already exists: {out.name}")
        return pd.read_csv(out)

    print(f"[4] ADMET filter (n={len(df_lip)})...")
    if not RDKIT_OK:
        df_lip.to_csv(out, index=False)
        return df_lip

    # PAINS filter
    pains_params = FilterCatalog.FilterCatalogParams()
    pains_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    pains_catalog = FilterCatalog.FilterCatalog(pains_params)

    rows = []
    for _, row in df_lip.iterrows():
        mol = Chem.MolFromSmiles(row["smiles"])
        if mol is None: continue

        is_pains = pains_catalog.HasMatch(mol)

        # LogS (ESOL approximation)
        mw   = row["mw_calc"]
        logp = row["logp"]
        hbd  = row["hbd"]
        rings = rdMolDescriptors.CalcNumRings(mol)
        rb   = row["rotbond"]
        logs = (0.16 - 0.63*logp - 0.0062*mw + 0.066*rb - 0.74*rings)

        # GI absorption (simple rule: TPSA≤140 + MW≤500 + not PAINS)
        gi_high = row["tpsa"] <= 140 and row["mw_calc"] <= 500

        if is_pains: continue
        if not gi_high: continue
        if logs < -6: continue

        rows.append({**row.to_dict(),
                     "rdkit_logs_esol": round(logs, 3),
                     "rdkit_gi_absorption": "High",
                     "rdkit_is_pains": is_pains})

    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"    Passed: {len(df)} → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 5: ChEMBL target prediction
# ════════════════════════════════════════════════════════════

def chembl_similarity(smiles, threshold=60):
    """Search ChEMBL for similar compounds and extract target hits."""
    r = requests.get(
        "https://www.ebi.ac.uk/chembl/api/data/similarity",
        params={"smiles": smiles, "similarity": threshold, "format": "json",
                "limit": 100},
        timeout=30
    )
    if r.status_code != 200: return []
    return r.json().get("molecules", {}).get("molecules", [])

def step5_targets(df_admet):
    out_pred = RES_DIR / "target_predictions.csv"
    out_cpd  = RES_DIR / "compounds_final.csv"
    if out_pred.exists() and out_cpd.exists():
        print(f"[5] Already exists: {out_pred.name}, {out_cpd.name}")
        return pd.read_csv(out_cpd), pd.read_csv(out_pred)

    print(f"[5] ChEMBL target prediction (n={len(df_admet)})...")
    pred_rows = []
    final_rows = []

    for i, row in df_admet.iterrows():
        molecules = chembl_similarity(row["smiles"], threshold=60)
        hits = []
        for mol in molecules:
            for activity in mol.get("activities", []):
                uniprot = activity.get("target_chembl_id")
                if uniprot not in TARGET_PROTEINS: continue
                pchembl = activity.get("pchembl_value")
                if pchembl is None or float(pchembl) < 5.0: continue
                hits.append({
                    "compound_name" : row["compound_name"],
                    "species"       : row.get("species",""),
                    "smiles"        : row["smiles"],
                    "cid"           : row.get("cid"),
                    "mw"            : row.get("mw_calc"),
                    "pchembl_value" : float(pchembl),
                    "uniprot_id"    : uniprot,
                    "protein_name"  : TARGET_PROTEINS[uniprot],
                })
        pred_rows.extend(hits)
        if hits:
            best = max(hits, key=lambda x: x["pchembl_value"])
            final_rows.append({**row.to_dict(),
                                "n_targets"  : len(set(h["uniprot_id"] for h in hits)),
                                "targets_hit": "; ".join(set(h["protein_name"] for h in hits)),
                                "best_pchembl": best["pchembl_value"]})
        if i % 20 == 0:
            print(f"    {i}/{len(df_admet)} processed, {len(final_rows)} with hits...")
        time.sleep(DELAY)

    df_pred  = pd.DataFrame(pred_rows)
    df_final = pd.DataFrame(final_rows)
    df_pred.to_csv(out_pred,  index=False)
    df_final.to_csv(out_cpd,  index=False)
    print(f"    Predictions: {len(df_pred)} rows → {out_pred.name}")
    print(f"    Final cpd  : {len(df_final)} → {out_cpd.name}")
    return df_final, df_pred

# ════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("="*55)
    print("Compound Filtering Pipeline")
    print("Tambunan et al. (2025)")
    print("="*55)
    df1 = step1_knapsack()
    df2 = step2_smiles(df1)
    df3 = step3_lipinski(df2)
    df4 = step4_admet(df3)
    df5, df_pred = step5_targets(df4)

    print("\n=== PIPELINE SUMMARY ===")
    print(f"  KNApSAcK raw        : {len(df1):>5}")
    print(f"  SMILES + dedup      : {len(df2):>5}")
    print(f"  Lipinski passed     : {len(df3):>5}")
    print(f"  ADMET passed        : {len(df4):>5}")
    print(f"  Final shortlist     : {len(df5):>5}")
    print(f"  Target predictions  : {len(df_pred):>5} rows")
