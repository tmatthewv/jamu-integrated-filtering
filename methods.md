# Detailed Methods Reference

This document supplements the Methods section of the published article with implementation details, parameter choices, and reproducibility notes.

---

## Stage 1 — Ethnopharmacology Data Collection

### KNApSAcK Jamu Database Scraping

**Script:** `scripts/01_ethnopharmacology/scrape_jamu.py`

The pipeline runs four sequential steps:

| Step | Function | Output |
|---|---|---|
| 1 | Parse product list from `top.php` | `jamu_products_raw.csv` |
| 2 | Scrape ingredients per product via `haigou.php` | `jamu_ingredients_raw.csv` |
| 3 | Resolve names via GBIF API | `jamu_ingredients_resolved.csv` |
| 4 | Aggregate frequencies | `jamu_frequency_final.csv` |

**GBIF resolution parameters:**
- Endpoint: `https://api.gbif.org/v1/species/match`
- Parameters: `kingdom=Plantae`, `strict=False`
- Fallback: secondary lookup via `usageKey` for SYNONYM matches
- Confidence threshold: ≥ 60 (below = UNRESOLVED)

**Taxonomy fixes** (`fix_unknown.py`):
| Raw name | Resolved to | Reason |
|---|---|---|
| Piperis Albi | Piper nigrum | Indonesian simplicia name |
| Atractylodis | Atractylodes macrocephala | Genus-only entry |
| Pinella | Pinellia ternata | Genus-only entry |

### Frequency Scoring

```
final_score = total_freq × weight

weight = 0.7 (Multi-claim plants — present in >1 category)
weight = 1.0 (Single-claim plants)

total_freq = MALE_FREQ + FEM_FREQ + NEU_FREQ
```

### Priority Plant Selection

```
Tukey upper fence = Q1 + 1.5 × IQR
Applied to final_score distribution (n=70 plants)
Threshold = 5.5
Result: 41 plants from 26 families (max 3 plants/family)
```

---

## Stage 2 — Visualization

### Correspondence Analysis (CA)

**Script:** `scripts/02_visualization/CA_tanaman.R`

- Package: FactoMineR (Lê et al., 2008)
- Input: contingency matrix 65 plants × 3 categories (total_freq ≥ 5)
- Dimensions retained: 2 (Dim1 = 61.3%, Dim2 = 38.7% of inertia)
- Label repulsion: ggrepel (force=3)

**CA plant parts:** `scripts/02_visualization/CA_bagian.R`
- Input: 9 standardized plant parts × 3 claim categories
- Cell values: number of producers per combination

### Network Visualization

**Script:** `scripts/02_visualization/network_tanaman.R`

- Package: igraph + ggraph
- Layout: stress-majorization (`create_layout(g, layout="stress")`)
- Seed: `set.seed(42)` for reproducibility
- Node size: min-max rescaled per level

---

## Stage 3 — Compound Filtering

**Script:** `scripts/03_compounds/pipeline.py`

### Filter cascade

```
KNApSAcK raw        : 3,699 entries
↓ PubChem SMILES + canonical dedup
SMILES valid        : 2,246 unique compounds
↓ Lipinski RO5 + Veber
  MW ≤ 500 Da
  LogP ≤ 5
  HBD ≤ 5
  HBA ≤ 10
  TPSA ≤ 140 Å²
  RotBond ≤ 10
Lipinski passed     : 1,464
↓ ADMET (RDKit)
  PAINS = False
  GI absorption = High
  LogS ≥ −6
  Lipinski violations ≤ 1
ADMET passed        : 1,292
↓ ChEMBL target prediction
  Tanimoto similarity ≥ 60%
  pChEMBL ≥ 5.0
  Target in {P10275, P03372, Q16665, P19838}
Final shortlist     : 31 compounds
```

---

## Stage 4 — Molecular Docking

### Protein Preparation

1. FASTA sequences retrieved from UniProt
2. AlphaFold2 structures predicted via AlphaFold Server (alphafoldserver.com)
3. pLDDT-based trimming: residues with pLDDT < 70 removed
4. ChimeraX preparation: remove solvent/ligands → addh → addcharge
5. Conversion to PDBQT via custom Python script

### Binding Site Determination

| Protein | Method | Reference structure | Native ligand |
|---|---|---|---|
| ANDR | Crystal superposition | PDB: 2AM9 | Testosterone (CID: 6013) |
| ESR1 | Crystal superposition | PDB: 1GWR | Estradiol (CID: 5757) |
| HIF1A | CB-Dock2 | — | 1,4-DPCA (CID: 459803) |
| NFKB1 | Crystal superposition | PDB: 8TQD | Parthenolide (CID: 6473881) |

### Validated Grid Box Parameters

| Protein | Center X | Center Y | Center Z | Box (Å) |
|---|---|---|---|---|
| ANDR | 26.78 | 2.36 | 4.61 | 22×22×22 |
| ESR1 | −12.83 | 11.25 | −8.94 | 22×22×22 |
| HIF1A | −3.00 | −11.00 | −11.00 | 20×20×20 |
| NFKB1 | −1.96 | 14.65 | −30.24 | 22×22×22 |

### Docking Parameters

```
exhaustiveness : 8
num_modes      : 9
energy_range   : 3 kcal/mol
grid_spacing   : 0.375 Å (default)
scoring        : Vina (empirical)
```

### Ligand Preparation (RDKit + Meeko)

```python
mol = Chem.MolFromSmiles(smiles)
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol, AllChem.ETKDGv3(), randomSeed=42)
AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
# → .sdf → meeko 0.5.0 → .pdbqt
```

---

## Stage 5 — Jamu Formulation Score

**Script:** `scripts/05_formulation/jamu_formulation_score.R`

### Target Relevance Matrix

| Protein | Male | Female | General | Justification |
|---|---|---|---|---|
| ANDR (P10275) | 1.00 | 0.10 | 0.40 | Androgen receptor: testosterone/DHT, muscle mass |
| ESR1 (P03372) | 0.10 | 1.00 | 0.40 | Estrogen receptor α: hormonal, bone density, endurance |
| HIF1A (Q16665) | 0.70 | 0.70 | 1.00 | Hypoxia adaptation, VO2max, aerobic capacity |
| NFKB1 (P19838) | 0.80 | 0.70 | 0.90 | Anti-inflammatory, DOMS recovery, immune modulation |

### JFS Formula

```
JFS(plant, cat) = TS(plant,cat) × BAS(plant,cat) × fs_norm(plant)

TS(plant, cat)  = FREQ_cat / total_freq
BAS(plant, cat) = Σ norm_pChEMBL(cpd) × TRS(target, cat)
fs_norm(plant)  = final_score / max(final_score)
norm_pChEMBL    = (pChEMBL - 4) / (max_pChEMBL - 4)

Plants without compound data: JFS × 0.3 penalty
```

### Statistical Tests

| Test | Purpose | Function |
|---|---|---|
| Wilcoxon rank-sum | Pairwise JFS distribution | `wilcox.test()` |
| Kruskal–Wallis | Overall group difference | `kruskal.test()` |
| Permutation test | Top-5 significance | Custom, n=10,000 |
| Bootstrap 95% CI | JFS uncertainty | `boot()`, R=2,000 |

---

## Downloading Protein Structures

AlphaFold2 predictions are not included in this repository. To reproduce:

1. Open https://alphafoldserver.com
2. Submit FASTA sequences from UniProt:
   - ANDR: https://www.uniprot.org/uniprot/P10275.fasta
   - ESR1: https://www.uniprot.org/uniprot/P03372.fasta
   - HIF1A: https://www.uniprot.org/uniprot/Q16665.fasta
   - NFKB1: https://www.uniprot.org/uniprot/P19838.fasta
3. Download `*_model_0.cif` files
4. Place in `data/raw/proteins/`
5. Run `scripts/04_docking/setup_docking.py` for conversion and prep
