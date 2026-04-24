#!/usr/bin/env python3
"""
scrape_jamu.py — KNApSAcK Jamu Database scraping pipeline
Tambunan et al. (2025) — Stage 1: Ethnopharmacology

Usage (from project root):
    conda activate docking
    python scripts/01_ethnopharmacology/scrape_jamu.py

Output files (data/processed/):
    jamu_products_raw.csv
    jamu_ingredients_raw.csv
    jamu_ingredients_resolved.csv
    taxonomy_cache.csv
    jamu_frequency_final.csv
"""

import os, re, time, json, requests, pandas as pd
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[2]
OUT_DIR    = ROOT / "data" / "processed"
OUT_DIR.mkdir(parents=True, exist_ok=True)

JAMU_BASE  = "https://www.knapsackfamily.com/jamu"
GBIF_URL   = "https://api.gbif.org/v1/species"
DELAY      = 0.5   # seconds between requests

# ── Category keywords ────────────────────────────────────────
MALE_KW    = ["lelaki","pria","jantan","perkasa","kejantanan","kelaki",
              "kaum adam","stamina pria","vitalitas pria"]
FEMALE_KW  = ["wanita","perempuan","kewanitaan","putri","feminim",
              "ibu","remaja putri","vitalitas wanita"]
GENERAL_KW = ["raga","badan","bugar","stamina","sehat","fitness",
              "kesehatan","kebugaran","prima","vitalitas"]

def classify_product(name: str) -> str:
    name_lower = name.lower()
    is_male    = any(k in name_lower for k in MALE_KW)
    is_female  = any(k in name_lower for k in FEMALE_KW)
    if is_male and not is_female:   return "Male-specific"
    if is_female and not is_male:   return "Female-specific"
    return "General health"

# ════════════════════════════════════════════════════════════
# STEP 1: Scrape product list
# ════════════════════════════════════════════════════════════

def step1_products():
    out = OUT_DIR / "jamu_products_raw.csv"
    if out.exists():
        print(f"[1] Already exists: {out.name}")
        return pd.read_csv(out)

    print("[1] Scraping product list (keyword: sehat)...")
    # Implementation: POST to top.php with keyword=sehat
    # Returns HTML table of product names
    products = []
    # [scraping logic omitted for brevity — see full implementation]
    df = pd.DataFrame(products, columns=["product_name","category"])
    df.to_csv(out, index=False)
    print(f"    Saved {len(df)} products → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 2: Scrape ingredients per product
# ════════════════════════════════════════════════════════════

def step2_ingredients(df_products):
    out = OUT_DIR / "jamu_ingredients_raw.csv"
    if out.exists():
        print(f"[2] Already exists: {out.name}")
        return pd.read_csv(out)

    print(f"[2] Scraping ingredients for {len(df_products)} products...")
    rows = []
    for i, row in df_products.iterrows():
        # POST to haigou.php with product name
        r = requests.post(f"{JAMU_BASE}/haigou.php",
                          data={"word": row["product_name"], "mode": "list"},
                          timeout=30)
        # Parse HTML table — extract: scientific name, vernacular, part, manufacturer
        # [parsing logic]
        time.sleep(DELAY)
        if i % 50 == 0:
            print(f"    {i}/{len(df_products)} products processed...")

    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"    Saved {len(df)} ingredient entries → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 3: GBIF taxonomy resolution
# ════════════════════════════════════════════════════════════

GBIF_CACHE = OUT_DIR / "taxonomy_cache.csv"

def resolve_gbif(name: str, cache: dict) -> dict:
    if name in cache:
        return cache[name]

    # Clean authorship citations
    clean = re.sub(r'\s+[A-Z][a-z]*\.?\s+.*$', '', name).strip()

    r = requests.get(f"{GBIF_URL}/match",
                     params={"name": clean, "kingdom": "Plantae", "strict": "false"},
                     timeout=15)
    if r.status_code != 200:
        return {"accepted_name": name, "status": "ERROR", "confidence": 0}

    data = r.json()
    if data.get("matchType") == "NONE":
        return {"accepted_name": name, "status": "UNRESOLVED", "confidence": 0}

    # If synonym, look up accepted name
    accepted = data.get("scientificName", name)
    if data.get("synonym", False) and data.get("usageKey"):
        r2 = requests.get(f"{GBIF_URL}/{data['usageKey']}", timeout=15)
        if r2.status_code == 200:
            accepted = r2.json().get("scientificName", accepted)

    result = {
        "input_name"   : name,
        "accepted_name": re.sub(r'\s+[A-Z].*$', '', accepted).strip(),
        "family"       : data.get("family", "Unknown"),
        "status"       : data.get("matchType", "UNKNOWN"),
        "confidence"   : data.get("confidence", 0)
    }
    cache[name] = result
    time.sleep(0.3)
    return result

def step3_taxonomy(df_ingredients):
    out = OUT_DIR / "jamu_ingredients_resolved.csv"
    if out.exists():
        print(f"[3] Already exists: {out.name}")
        return pd.read_csv(out)

    print("[3] Resolving taxonomy via GBIF...")
    cache = {}
    if GBIF_CACHE.exists():
        cache_df = pd.read_csv(GBIF_CACHE)
        cache = cache_df.set_index("input_name").to_dict(orient="index")

    unique_names = df_ingredients["scientific_name"].dropna().unique()
    print(f"    Querying {len(unique_names)} unique names...")

    results = {}
    for i, name in enumerate(unique_names):
        results[name] = resolve_gbif(name, cache)
        if i % 20 == 0:
            print(f"    {i}/{len(unique_names)} resolved...")
            pd.DataFrame(list(results.values())).to_csv(GBIF_CACHE, index=False)

    pd.DataFrame(list(results.values())).to_csv(GBIF_CACHE, index=False)

    df = df_ingredients.copy()
    df["accepted_name"] = df["scientific_name"].map(
        lambda x: results.get(x, {}).get("accepted_name", x))
    df["family"] = df["scientific_name"].map(
        lambda x: results.get(x, {}).get("family", "Unknown"))
    df.to_csv(out, index=False)
    print(f"    Saved {len(df)} resolved entries → {out.name}")
    return df

# ════════════════════════════════════════════════════════════
# STEP 4: Frequency aggregation + scoring
# ════════════════════════════════════════════════════════════

def step4_frequency(df_resolved):
    out = OUT_DIR / "jamu_frequency_final.csv"
    if out.exists():
        print(f"[4] Already exists: {out.name}")
        return pd.read_csv(out)

    print("[4] Aggregating frequency scores...")

    freq = (df_resolved
            .groupby(["accepted_name","family","category"])["manufacturer"]
            .nunique()
            .reset_index(name="freq"))

    pivot = freq.pivot_table(
        index=["accepted_name","family"],
        columns="category",
        values="freq",
        fill_value=0
    ).reset_index()

    for col in ["Male-specific","Female-specific","General health"]:
        if col not in pivot.columns:
            pivot[col] = 0

    pivot.columns.name = None
    pivot = pivot.rename(columns={
        "accepted_name"  : "TANAMAN",
        "Male-specific"  : "MALE_FREQ",
        "Female-specific": "FEM_FREQ",
        "General health" : "NEU_FREQ"
    })

    pivot["total_freq"] = pivot["MALE_FREQ"] + pivot["FEM_FREQ"] + pivot["NEU_FREQ"]

    def get_group(row):
        cats = sum([row["MALE_FREQ"]>0, row["FEM_FREQ"]>0, row["NEU_FREQ"]>0])
        if cats > 1: return "Mixed"
        if row["MALE_FREQ"] > 0: return "Male"
        if row["FEM_FREQ"]  > 0: return "Female"
        return "Neutral"

    pivot["group"]       = pivot.apply(get_group, axis=1)
    pivot["weight"]      = pivot["group"].map({"Mixed":0.7}).fillna(1.0)
    pivot["final_score"] = pivot["total_freq"] * pivot["weight"]

    pivot.to_csv(out, index=False)
    print(f"    Saved {len(pivot)} plants → {out.name}")
    return pivot

# ════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("="*55)
    print("KNApSAcK Jamu Scraping Pipeline")
    print("Tambunan et al. (2025)")
    print("="*55)
    df_products    = step1_products()
    df_ingredients = step2_ingredients(df_products)
    df_resolved    = step3_taxonomy(df_ingredients)
    df_frequency   = step4_frequency(df_resolved)
    print("\nPipeline selesai.")
    print(f"  Products  : {len(df_products)}")
    print(f"  Ingredients: {len(df_ingredients)}")
    print(f"  Species   : {df_resolved['accepted_name'].nunique()}")
    print(f"  Final freq: {len(df_frequency)} plants")
