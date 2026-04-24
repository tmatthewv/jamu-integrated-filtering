# Protein Target Selection Rationale

This document provides the biological justification for the four protein targets selected for molecular docking and Jamu Formulation Score calculation.

---

## Selection Workflow

Targets were identified through three complementary approaches:

1. **GeneCards queries** on: "immune modulation", "androgen signaling", "estrogen signaling", "ATP production", "muscle endurance"
2. **Manual inclusion** of sex-hormone receptors (ANDR, ESR1) based on physiological relevance to sex-specific stamina
3. **Reactome pathway enrichment** via STRING to validate mechanistic relevance

Six candidate targets were initially identified. IL6RA (P08887) and IGF1R (P08069) were excluded due to no matching compounds in the filtered dataset.

---

## Selected Targets

### 1. Androgen Receptor — ANDR (P10275)

**Relevance to stamina:** The androgen receptor mediates testosterone and dihydrotestosterone signaling, which directly regulates skeletal muscle mass, erythropoiesis, and exercise-induced anabolic responses. Testosterone acts through ANDR to promote nitrogen retention and muscle protein synthesis, both critical for sustained physical performance in males.

**Rationale for inclusion:** Manually included. Associated with male-specific stamina maintenance and matched to ≥1 compound in the dataset (Estriol, Eugenol, 3-Hydroxyflavone).

**Key pathways:** Androgen receptor signaling, nuclear receptor transcription, SUMOylation of intracellular receptors.

**TRS:** Male = 1.00, Female = 0.10, General = 0.40

---

### 2. Estrogen Receptor Alpha — ESR1 (P03372)

**Relevance to stamina:** ESR1 mediates the effects of 17β-estradiol on skeletal muscle, bone mineral density, and cardiovascular adaptation. Estrogen signaling through ESR1 modulates mitochondrial biogenesis, lipid metabolism during endurance exercise, and post-exercise inflammatory resolution — effects that are particularly prominent in females.

**Rationale for inclusion:** Manually included. Associated with female-specific stamina and matched to multiple compounds (Oestrone, Equol, Estriol, Chrysin, Daidzein, Genistein, Glabranin, Apigenin, Butylparaben, alpha-Zearalenol, Isoglabranin).

**Key pathways:** Nuclear receptor transcription, estrogen receptor signaling, mitochondrial biogenesis.

**TRS:** Male = 0.10, Female = 1.00, General = 0.40

---

### 3. Hypoxia-Inducible Factor 1-alpha — HIF1A (Q16665)

**Relevance to stamina:** HIF-1α is the master regulator of cellular responses to hypoxia and is centrally involved in aerobic capacity, VO2max, and endurance performance. It transcriptionally activates genes involved in erythropoiesis (EPO), glycolysis, angiogenesis, and mitochondrial function. HIF-1α activation during exercise is considered a universal mechanism for adaptation across all athletes.

**Rationale for inclusion:** Identified through GeneCards ("ATP production", "muscle endurance"). Matched to 8 compounds (Estradiol, Hesperetin, L-Tryptophan, Piperlongumine, Furaltadone, Adenosine, L-Glutamine, Se-Methylselenocysteine*).

**Key pathways:** Hypoxia response, EPO signaling, glycolysis/gluconeogenesis, VEGF signaling.

**TRS:** Male = 0.70, Female = 0.70, General = 1.00

*Se-Methylselenocysteine excluded from docking analysis.

---

### 4. Nuclear Factor NF-κB p105 — NFKB1 (P19838)

**Relevance to stamina:** NF-κB1 is a transcription factor that controls the inflammatory response to exercise-induced muscle damage. It regulates the balance between glycolysis and mitochondrial oxidative phosphorylation (OXPHOS), and its modulation is relevant to reducing delayed-onset muscle soreness (DOMS) and accelerating recovery. NF-κB1 inhibition has been associated with improved exercise tolerance across sexes.

**Rationale for inclusion:** Identified through GeneCards ("immune modulation"). Matched to 11 compounds (Withanolide, Gibberellin A3, Homobutein, Isoliquiritigenin, Biochanin A, Prunetin, 2-Hydroxycinnamic acid, Tryptamine, Agmatine, 5-Aminopentanoic acid, Se-Methylselenocysteine*).

**Key pathways:** NF-κB signaling, inflammatory response, OXPHOS regulation, deubiquitination.

**TRS:** Male = 0.80, Female = 0.70, General = 0.90

---

## Excluded Targets

### IL-6 Receptor — IL6RA (P08887)
Identified through GeneCards ("immune modulation"). Excluded because no compounds in the filtered dataset matched this target at the required Tanimoto similarity and pChEMBL threshold.

### IGF-1 Receptor — IGF1R (P08069)
Identified through GeneCards ("muscle endurance"). Excluded for the same reason as IL6RA — zero compound matches in the curated dataset.

---

## Pathway Enrichment Validation

Reactome pathway enrichment was performed via STRING (https://string-db.org) using the four UniProt accessions as input. Enriched pathways confirmed the mechanistic relevance of all four targets to athletic performance:

- SUMOylation of intracellular receptors (ANDR, ESR1)
- Nuclear receptor transcription pathways (ANDR, ESR1)
- Cellular responses to stress (HIF1A, NFKB1)
- Deubiquitination processes (NFKB1)
- Transcriptional regulation by VENTX (NFKB1)
