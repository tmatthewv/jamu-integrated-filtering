# Integrative Ethnopharmacological and In-Silico Strategy for Identifying Indonesian Jamu to Enhance Athletic Stamina

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![R 4.5+](https://img.shields.io/badge/R-4.5+-276DC3.svg)](https://www.r-project.org/)

---

## Overview

This repository contains all scripts, data, and reproducible workflows for an integrative study that combines ethnopharmacological analysis of Indonesian jamu formulations with in-silico drug discovery methods to identify plant-derived candidates for athletic stamina enhancement.

The pipeline proceeds through five stages:

```
KNApSAcK Jamu Database (5,310 products)
        ↓  [01_ethnopharmacology]
41 priority plants · 26 families
        ↓  [02_visualization]
Network · CA · PRISMA funnel
        ↓  [03_compounds]
31 candidate compounds (drug-likeness + ADMET + ChEMBL)
        ↓  [04_docking]
Molecular docking → 4 targets (ANDR, ESR1, HIF1A, NFKB1)
        ↓  [05_formulation]
Jamu Formulation Score (JFS) → Male / Female / General
```

---

## Repository Structure

```
jamu-stamina-insilico/
├── README.md
├── LICENSE
├── environment.yml              # conda environment (Python)
├── renv.lock                    # R package lockfile
├── .gitignore
│
├── data/
│   ├── raw/                     # Input data (read-only)
│   │   ├── bobot_tanaman.csv    # 70 plants with frequency scores
│   │   ├── JAMU1.csv            # Original jamu dataset (JAMU1)
│   │   └── target_relevance_matrix.csv
│   ├── processed/               # Intermediate outputs
│   │   ├── jamu_products_raw.csv
│   │   ├── jamu_ingredients_raw.csv
│   │   ├── jamu_ingredients_resolved.csv
│   │   ├── taxonomy_cache.csv
│   │   ├── jamu_frequency_final.csv
│   │   ├── knapsack_raw.csv
│   │   ├── compounds_with_smiles.csv
│   │   ├── compounds_lipinski.csv
│   │   └── compounds_admet.csv
│   └── results/                 # Final outputs
│       ├── compounds_final.csv
│       ├── target_predictions.csv
│       ├── docking_results_all.csv
│       ├── docking_results_best.csv
│       ├── jamu_formulation_scores.csv
│       └── jamu_formulation_recommended.csv
│
├── scripts/
│   ├── 01_ethnopharmacology/
│   │   ├── scrape_jamu.py       # KNApSAcK scraping + GBIF resolution
│   │   └── fix_unknown.py       # Manual taxonomy fixes
│   ├── 02_visualization/
│   │   ├── network_tanaman.R    # Plant-family network
│   │   ├── CA_tanaman.R         # Correspondence analysis (species)
│   │   ├── CA_bagian.R          # Correspondence analysis (plant parts)
│   │   ├── network_top.R        # Top-35 radial network
│   │   ├── network_compound_plant_family.R  # Compound-plant-family network
│   │   ├── prisma_funnel.R      # PRISMA filtering diagram
│   │   └── funnel_report.R      # Filter report (Lipinski + ADMET detail)
│   ├── 03_compounds/
│   │   └── pipeline.py          # Full compound filtering pipeline
│   ├── 04_docking/
│   │   ├── setup_docking.py     # Grid box config generation
│   │   ├── run_docking.sh       # Batch AutoDock Vina runner
│   │   └── extract_results.py  # Parse docking logs → CSV
│   └── 05_formulation/
│       └── jamu_formulation_score.R  # JFS calculation + statistics
│
├── docs/
│   ├── methods.md               # Detailed methods reference
│   └── target_justification.md # Protein target selection rationale
│
└── figures/                     # Publication-ready figures (PNG, 300 dpi)
    ├── fig1_prisma_funnel.png
    ├── fig2_venn_targets.png
    ├── fig3_reactome_pathway.png
    ├── fig4_CA_bagian.png
    ├── fig5_CA_tanaman.png
    ├── fig6_network_compound_plant_family.png
    ├── supp_network_tanaman.png
    ├── supp_funnel_report.png
    └── supp_jfs_heatmap.png
```

---

## Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/tmatthewv/jamu-stamina-insilico.git
cd jamu-stamina-insilico
```

### 2. Set up Python environment

```bash
conda env create -f environment.yml
conda activate docking
```

### 3. Set up R environment

```r
install.packages("renv")
renv::restore()
```

### 4. Run the pipeline

Each stage can be run independently. Run in order for full reproducibility:

```bash
# Stage 1: Ethnopharmacology
cd scripts/01_ethnopharmacology
python scrape_jamu.py

# Stage 2: Visualization (in RStudio, set working dir to project root)
# Rscript scripts/02_visualization/CA_tanaman.R
# Rscript scripts/02_visualization/network_tanaman.R
# Rscript scripts/02_visualization/prisma_funnel.R

# Stage 3: Compound pipeline
cd scripts/03_compounds
python pipeline.py

# Stage 4: Molecular docking
cd scripts/04_docking
python setup_docking.py
bash run_docking.sh
python extract_results.py

# Stage 5: Formulation scoring
# Rscript scripts/05_formulation/jamu_formulation_score.R
```

---

## Dependencies

### Python (≥ 3.10)
| Package | Version | Purpose |
|---|---|---|
| rdkit | ≥ 2023.09 | Molecular descriptors, 3D conformers |
| pandas | ≥ 2.0 | Data manipulation |
| requests | ≥ 2.31 | API calls (PubChem, GBIF, ChEMBL) |
| meeko | 0.5.0 | PDBQT preparation |
| biopython | ≥ 1.81 | Sequence handling |

### R (≥ 4.5)
| Package | Purpose |
|---|---|
| FactoMineR | Correspondence analysis |
| ggraph, igraph | Network visualization |
| ggplot2, ggrepel | Plotting |
| dplyr, tidyr | Data wrangling |
| boot | Bootstrap confidence intervals |
| patchwork | Multi-panel figures |

### External tools
| Tool | Version | Purpose |
|---|---|---|
| AutoDock Vina | 1.2.5 | Molecular docking |
| UCSF ChimeraX | 1.11 | Protein preparation, visualization |
| Open Babel | 3.1.0 | File format conversion |
| AlphaFold Server | — | Protein structure prediction |

---

## Methods Summary

### Ethnopharmacology (Stage 1)
- 5,310 KNApSAcK Jamu products scraped; 368 filtered by keyword "sehat"
- Classified into Male-specific (n=103), Female-specific (n=63), General health (n=202)
- 15,422 ingredient entries resolved to 234 unique accepted species via GBIF API
- Priority plants selected using Tukey upper fence (Q1 + 1.5×IQR = 5.5) with max 3 plants/family → **41 plants, 26 families**

### Compound Filtering (Stage 3)
- 3,699 raw KNApSAcK compounds → 2,246 (SMILES + dedup) → 1,464 (Lipinski/Veber) → 1,292 (ADMET) → **31 final** (ChEMBL target prediction, pChEMBL ≥ 5.0)

### Molecular Docking (Stage 4)
- Protein structures: AlphaFold2 predictions, pLDDT-trimmed (< 70 removed)
- Binding sites: co-crystallized ligand superposition (ANDR, ESR1, NFKB1) + CB-Dock2 (HIF1A)
- AutoDock Vina v1.2.5, exhaustiveness=8, 9 poses, energy range=3 kcal/mol

### Jamu Formulation Score (Stage 5)
```
JFS(plant, category) = TS × BAS × fs_norm

TS  = tendency score (FREQ_cat / total_freq)
BAS = Σ norm_pChEMBL × target_relevance_score
fs_norm = final_score / max(final_score)
```
- Statistical validation: Wilcoxon rank-sum, Kruskal–Wallis, permutation test (n=10,000), bootstrap 95% CI

---

## Data Availability

Raw scraping outputs and intermediate files are included in `data/processed/`. Protein structure files (AlphaFold2 `.cif` and `.pdb`) are not included due to size; see `docs/methods.md` for download instructions.

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

---

## Contact

**Matthew Valentino Tambunan** — tambunan.matthewv@gmail.com  
Department of Pharmacy, Faculty of Mathematics and Natural Sciences, Udayana University, Bali, Indonesia
