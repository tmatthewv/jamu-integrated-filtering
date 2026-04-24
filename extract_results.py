# ============================================================
# Extract docking results → summary CSV + ranking
# ============================================================

import os, re, pandas as pd

RESULT_DIR = "results"
PROTEINS   = ["2AM9", "3ERT", "3HQR", "1SVC"]
PROT_NAMES = {
    "2AM9": "Androgen Receptor (AR)",
    "3ERT": "Estrogen Receptor alpha (ESR1)",
    "3HQR": "HIF-1alpha (HIF1A)",
    "1SVC": "NF-kB1 (NFKB1)",
}

rows = []
for prot in PROTEINS:
    prot_dir = os.path.join(RESULT_DIR, prot)
    if not os.path.isdir(prot_dir):
        continue
    for logfile in os.listdir(prot_dir):
        if not logfile.endswith("_log.txt"):
            continue
        name = logfile.replace("_log.txt", "")
        logpath = os.path.join(prot_dir, logfile)
        with open(logpath) as f:
            content = f.read()
        # Ambil semua mode
        modes = re.findall(r'^\s+(\d+)\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)',
                           content, re.MULTILINE)
        if modes:
            best = modes[0]
            rows.append({
                "compound"    : name,
                "protein"     : PROT_NAMES.get(prot, prot),
                "pdb_id"      : prot,
                "mode_1_affinity": float(best[1]),
                "mode_1_rmsd_lb" : float(best[2]),
                "mode_1_rmsd_ub" : float(best[3]),
                "n_modes"        : len(modes),
            })
        else:
            print(f"  [WARN] No modes parsed: {logfile}")

if not rows:
    print("Tidak ada hasil — pastikan docking sudah selesai")
else:
    df = pd.DataFrame(rows)

    # Merge dengan compound info
    if os.path.exists("compounds_final.csv"):
        df_cpd = pd.read_csv("compounds_final.csv")[
            ["compound_name","smiles","mw_calc","logp","found_in_species",
             "n_species","best_pchembl","targets_hit"]
        ]
        df_cpd["compound"] = df_cpd["compound_name"].str.replace(
            r'[^a-zA-Z0-9_-]', '_', regex=True).str[:40]
        df = df.merge(df_cpd, on="compound", how="left")

    df = df.sort_values("mode_1_affinity")  # lebih negatif = lebih baik

    df.to_csv("docking_results_all.csv", index=False)
    print(f"Disimpan: docking_results_all.csv ({len(df)} baris)")

    # Top 10 per protein
    print("\n=== TOP 10 PER PROTEIN (binding affinity, kcal/mol) ===")
    for prot_name in df["protein"].unique():
        sub = df[df["protein"]==prot_name].head(10)
        print(f"\n{prot_name}:")
        print(sub[["compound","mode_1_affinity","best_pchembl",
                    "found_in_species"]].to_string(index=False))

    # Summary terbaik per senyawa (lintas protein)
    best_per_cpd = df.loc[df.groupby("compound")["mode_1_affinity"].idxmin()]
    best_per_cpd = best_per_cpd.sort_values("mode_1_affinity")
    best_per_cpd.to_csv("docking_results_best.csv", index=False)
    print(f"\nTop 10 senyawa lintas semua protein:")
    print(best_per_cpd[["compound","protein","mode_1_affinity",
                          "best_pchembl","mw_calc"]].head(10).to_string(index=False))
