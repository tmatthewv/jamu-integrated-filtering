#!/usr/bin/env python3
"""
fix_unknown.py — Manual taxonomy fixes for non-standard simplicia names
Tambunan et al. (2025)

Some entries in KNApSAcK use Indonesian simplicia names (pharmaceutical
Latin-like names) rather than standard binomial nomenclature.
These cannot be resolved by GBIF and require manual mapping.

Usage:
    python scripts/01_ethnopharmacology/fix_unknown.py
"""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PRO  = ROOT / "data" / "processed"

MANUAL_FIXES = {
    "Piperis Albi"              : ("Piper nigrum",             "Piperaceae"),
    "Atractylodis"              : ("Atractylodes macrocephala","Asteraceae"),
    "Pinella"                   : ("Pinellia ternata",         "Araceae"),
    "Zingiberis Officinalis"    : ("Zingiber officinale",      "Zingiberaceae"),
    "Curcumae Xanthorrhizae"    : ("Curcuma xanthorrhiza",     "Zingiberaceae"),
    "Languatis"                 : ("Alpinia galanga",          "Zingiberaceae"),
    "Kaempferiae Galangae"      : ("Kaempferia galanga",       "Zingiberaceae"),
    "Caryophylli"               : ("Syzygium aromaticum",      "Myrtaceae"),
}

def apply_fixes():
    f = PRO / "jamu_ingredients_resolved.csv"
    if not f.exists():
        print("ERROR: Run scrape_jamu.py first")
        return

    df = pd.read_csv(f)
    fixes_applied = 0

    for raw, (accepted, family) in MANUAL_FIXES.items():
        mask = df["scientific_name"].str.contains(raw, case=False, na=False)
        if mask.sum() > 0:
            df.loc[mask, "accepted_name"] = accepted
            df.loc[mask, "family"]        = family
            fixes_applied += mask.sum()
            print(f"  Fixed: '{raw}' → '{accepted}' ({mask.sum()} rows)")

    df.to_csv(f, index=False)
    print(f"\nApplied {fixes_applied} manual taxonomy fixes.")

if __name__ == "__main__":
    apply_fixes()
