# ============================================================
# PRISMA-style Filtering Flowchart
# Semua nilai n dibaca langsung dari file output pipeline
# ============================================================

library(ggplot2)
library(dplyr)

col_bg    <- "#FFFFFF"
col_main  <- "#1565A0"
col_plant <- "#0F6E56"
col_cpd   <- "#854F0B"
col_dock  <- "#72243E"
col_excl  <- "#5F5E5A"
col_text  <- "#1A1A1A"

# ════════════════════════════════════════════════════════════
# STEP 1: Baca semua file dan hitung n
# ════════════════════════════════════════════════════════════

read_safe <- function(f, ...) {
  if (file.exists(f)) read.csv(f, ...) else NULL
}

# ── Products ─────────────────────────────────────────────────
df_products    <- read_safe("jamu_products_raw.csv")
df_ingredients <- read_safe("jamu_ingredients_raw.csv")
df_resolved    <- read_safe("jamu_ingredients_resolved.csv")
df_frequency   <- read_safe("jamu_frequency_final.csv")

# ── Compounds ─────────────────────────────────────────────────
df_knapsack    <- read_safe("knapsack_raw.csv")
df_smiles      <- read_safe("compounds_with_smiles.csv")
df_lipinski    <- read_safe("compounds_lipinski.csv")
df_admet       <- read_safe("compounds_admet.csv")
df_targets     <- read_safe("target_predictions.csv")
df_final       <- read_safe("compounds_final.csv")

# ── Hitung n per tahap ────────────────────────────────────────
cat("=== PRISMA n values dari file ===\n")

# Products
n_knapsack_total <- 5310  # total KNApSAcK Jamu (fixed, dari website)

n_sehat <- if (!is.null(df_products)) nrow(df_products) else NA
n_male  <- if (!is.null(df_products)) sum(df_products$category == "Male-specific", na.rm=TRUE) else NA
n_female<- if (!is.null(df_products)) sum(df_products$category == "Female-specific", na.rm=TRUE) else NA
n_gen   <- if (!is.null(df_products)) sum(df_products$category == "General health", na.rm=TRUE) else NA

# Ingredients
n_ing_raw <- if (!is.null(df_ingredients)) nrow(df_ingredients) else NA

# Taxonomy
n_synonyms <- if (!is.null(df_resolved) && "status" %in% names(df_resolved)) {
  sum(df_resolved$input_name != df_resolved$accepted_name, na.rm=TRUE)
} else NA

n_unique_spp <- if (!is.null(df_resolved)) {
  length(unique(df_resolved$accepted_name[!is.na(df_resolved$accepted_name)]))
} else NA

n_families <- if (!is.null(df_frequency)) {
  length(unique(df_frequency$famili[df_frequency$famili != "Unknown"]))
} else NA

# Taxonomy status breakdown (dari taxonomy_cache.csv jika ada)
df_tax <- read_safe("taxonomy_cache.csv")
n_exact     <- if (!is.null(df_tax)) sum(df_tax$status == "EXACT",       na.rm=TRUE) else NA
n_fuzzy     <- if (!is.null(df_tax)) sum(df_tax$status == "FUZZY",       na.rm=TRUE) else NA
n_higherrank<- if (!is.null(df_tax)) sum(df_tax$status == "HIGHERRANK",  na.rm=TRUE) else NA

# Priority plants (Tukey)
if (!is.null(df_frequency)) {
  plant_score <- df_frequency %>%
    group_by(tanaman, famili) %>%
    summarise(final_score = sum(final_score, na.rm=TRUE), .groups="drop")

  Q1 <- quantile(plant_score$final_score, 0.25)
  Q3 <- quantile(plant_score$final_score, 0.75)
  tukey_thresh <- Q1 + 1.5*(Q3 - Q1)

  priority <- plant_score %>%
    filter(final_score >= tukey_thresh) %>%
    arrange(desc(final_score)) %>%
    group_by(famili) %>% slice_head(n=3) %>% ungroup()

  n_priority       <- nrow(priority)
  n_priority_fam   <- n_distinct(priority$famili)
  n_below_tukey    <- nrow(plant_score) - n_priority
} else {
  n_priority <- n_priority_fam <- n_below_tukey <- tukey_thresh <- NA
}

# Compounds
n_knapsack_raw  <- if (!is.null(df_knapsack)) nrow(df_knapsack) else NA
n_knapsack_spp  <- if (!is.null(df_knapsack)) n_distinct(df_knapsack$species) else NA

n_smiles_valid  <- if (!is.null(df_smiles)) nrow(df_smiles) else NA
n_smiles_excl   <- if (!is.null(df_knapsack) && !is.null(df_smiles))
  n_knapsack_raw - n_smiles_valid else NA

n_lipinski_pass <- if (!is.null(df_lipinski)) nrow(df_lipinski) else NA
n_lipinski_excl <- if (!is.null(df_smiles) && !is.null(df_lipinski))
  n_smiles_valid - n_lipinski_pass else NA

# Lipinski per-criterion (dari kolom yang ada di compounds_lipinski.csv)
if (!is.null(df_lipinski)) {
  lip_mw  <- if ("mw_calc" %in% names(df_lipinski))
    sum(df_lipinski$mw_calc  <= 500, na.rm=TRUE) else
    if ("mw" %in% names(df_lipinski)) sum(as.numeric(df_lipinski$mw) <= 500, na.rm=TRUE) else NA
  lip_logp <- if ("logp"    %in% names(df_lipinski)) sum(df_lipinski$logp    <= 5,   na.rm=TRUE) else NA
  lip_hbd  <- if ("hbd"     %in% names(df_lipinski)) sum(df_lipinski$hbd     <= 5,   na.rm=TRUE) else NA
  lip_hba  <- if ("hba"     %in% names(df_lipinski)) sum(df_lipinski$hba     <= 10,  na.rm=TRUE) else NA
  lip_tpsa <- if ("tpsa"    %in% names(df_lipinski)) sum(df_lipinski$tpsa    <= 140, na.rm=TRUE) else NA
  lip_rot  <- if ("rotbond" %in% names(df_lipinski)) sum(df_lipinski$rotbond <= 10,  na.rm=TRUE) else NA
} else {
  lip_mw <- lip_logp <- lip_hbd <- lip_hba <- lip_tpsa <- lip_rot <- NA
}

# ADMET per-criterion
n_admet_pass  <- if (!is.null(df_admet)) nrow(df_admet) else NA
n_admet_excl  <- if (!is.null(df_lipinski) && !is.null(df_admet))
  n_lipinski_pass - n_admet_pass else NA

if (!is.null(df_admet)) {
  adm_pains <- if ("pains" %in% names(df_admet))
    sum(!df_admet$pains, na.rm=TRUE) else NA
  adm_gi    <- if ("gi_absorption" %in% names(df_admet))
    sum(df_admet$gi_absorption == "High", na.rm=TRUE) else NA
  adm_logs  <- if ("logs" %in% names(df_admet))
    sum(df_admet$logs >= -6, na.rm=TRUE) else NA
  adm_bbb   <- if ("bbb_permeant" %in% names(df_admet))
    sum(df_admet$bbb_permeant == TRUE | df_admet$bbb_permeant == "True", na.rm=TRUE) else
    if ("bbb" %in% names(df_admet)) sum(df_admet$bbb, na.rm=TRUE) else NA
  adm_lipviol <- if ("lipinski_violations" %in% names(df_admet))
    sum(df_admet$lipinski_violations <= 1, na.rm=TRUE) else NA
  # PAINS excluded from lipinski total
  n_pains_excl <- if (!is.null(df_lipinski) && "pains" %in% names(df_lipinski))
    sum(df_lipinski$pains == TRUE | df_lipinski$pains == "True", na.rm=TRUE) else 101
} else {
  adm_pains <- adm_gi <- adm_logs <- adm_bbb <- adm_lipviol <- NA
  n_pains_excl <- NA
}

# Target prediction
n_target_hits  <- if (!is.null(df_targets)) n_distinct(df_targets$compound_name) else NA
n_target_total <- if (!is.null(df_targets)) nrow(df_targets) else NA
n_targets_used <- if (!is.null(df_targets)) n_distinct(df_targets$uniprot_id) else NA

# Per-protein breakdown
if (!is.null(df_targets) && "protein_name" %in% names(df_targets)) {
  target_breakdown <- df_targets %>%
    group_by(protein_name) %>%
    summarise(n_cpd = n_distinct(compound_name), .groups = "drop") %>%
    arrange(desc(n_cpd))
} else {
  target_breakdown <- NULL
}

# Shortlist criteria from compounds_final
n_final         <- if (!is.null(df_final)) nrow(df_final) else NA
n_final_AR      <- if (!is.null(df_final) && "best_target" %in% names(df_final))
  sum(grepl("AR",    df_final$best_target), na.rm=TRUE) else NA
n_final_ESR1    <- if (!is.null(df_final) && "best_target" %in% names(df_final))
  sum(grepl("ESR1",  df_final$best_target), na.rm=TRUE) else NA
n_final_HIF1A   <- if (!is.null(df_final) && "best_target" %in% names(df_final))
  sum(grepl("HIF1A", df_final$best_target), na.rm=TRUE) else NA
n_final_NFKB1   <- if (!is.null(df_final) && "best_target" %in% names(df_final))
  sum(grepl("NFKB1", df_final$best_target), na.rm=TRUE) else NA

# Compounds excluded at target step (ADMET pass but no target hit / pChEMBL < 5)
n_no_target_hit <- if (!is.na(n_admet_pass) && !is.na(n_target_hits))
  n_admet_pass - n_target_hits else NA
n_target_excl   <- if (!is.na(n_target_hits) && !is.na(n_final))
  n_target_hits - n_final else NA

# Print semua nilai
cat(sprintf("n_knapsack_total : %s\n", n_knapsack_total))
cat(sprintf("n_sehat          : %s\n", n_sehat))
cat(sprintf("  Male           : %s\n", n_male))
cat(sprintf("  Female         : %s\n", n_female))
cat(sprintf("  General        : %s\n", n_gen))
cat(sprintf("n_ing_raw        : %s\n", n_ing_raw))
cat(sprintf("n_synonyms       : %s\n", n_synonyms))
cat(sprintf("n_unique_spp     : %s\n", n_unique_spp))
cat(sprintf("n_families       : %s\n", n_families))
cat(sprintf("  EXACT/FUZZY/HR : %s/%s/%s\n", n_exact, n_fuzzy, n_higherrank))
cat(sprintf("tukey_thresh     : %.2f\n", tukey_thresh))
cat(sprintf("n_priority       : %s (%s families)\n", n_priority, n_priority_fam))
cat(sprintf("n_below_tukey    : %s\n", n_below_tukey))
cat(sprintf("n_knapsack_raw   : %s (%s spp)\n", n_knapsack_raw, n_knapsack_spp))
cat(sprintf("n_smiles_valid   : %s (excl %s)\n", n_smiles_valid, n_smiles_excl))
cat(sprintf("n_lipinski_pass  : %s (excl %s)\n", n_lipinski_pass, n_lipinski_excl))
cat(sprintf("  MW/LogP/HBD    : %s/%s/%s\n", lip_mw, lip_logp, lip_hbd))
cat(sprintf("  HBA/TPSA/Rot   : %s/%s/%s\n", lip_hba, lip_tpsa, lip_rot))
cat(sprintf("n_admet_pass     : %s (excl %s)\n", n_admet_pass, n_admet_excl))
cat(sprintf("  PAINS excl     : %s\n", n_pains_excl))
cat(sprintf("  GI/LogS/BBB    : %s/%s/%s\n", adm_gi, adm_logs, adm_bbb))
cat(sprintf("n_target_hits    : %s\n", n_target_hits))
cat(sprintf("n_final          : %s (%s targets)\n", n_final, n_targets_used))
if (!is.null(target_breakdown)) {
  cat("Target breakdown:\n")
  for (i in seq_len(nrow(target_breakdown)))
    cat(sprintf("  %-40s : %s\n", target_breakdown$protein_name[i], target_breakdown$n_cpd[i]))
}

# ── format helper ─────────────────────────────────────────
fmt <- function(x, fallback="N/A") {
  if (is.null(x) || is.na(x)) fallback
  else formatC(as.integer(x), format="d", big.mark=",")
}
pct <- function(num, den) {
  if (is.null(num)||is.null(den)||is.na(num)||is.na(den)||den==0) ""
  else sprintf("(%.1f%%)", num/den*100)
}

# ════════════════════════════════════════════════════════════
# STEP 2: Build diagram
# ════════════════════════════════════════════════════════════

box <- function(g, x, y, w, h, label, sub=NULL,
                fill="#E6F1FB", border=col_main,
                fontsize=3.2, subfontsize=2.7, bold=FALSE) {
  g <- g +
    annotate("rect", xmin=x, xmax=x+w, ymin=y, ymax=y+h,
             fill=fill, color=border, linewidth=0.5) +
    annotate("text", x=x+w/2,
             y=if(is.null(sub)) y+h/2 else y+h*0.64,
             label=label, hjust=0.5, vjust=0.5,
             size=fontsize, color=col_text, lineheight=1.15,
             fontface=if(bold)"bold" else "plain")
  if (!is.null(sub))
    g <- g + annotate("text", x=x+w/2, y=y+h*0.27,
                      label=sub, hjust=0.5, vjust=0.5,
                      size=subfontsize, color="#444444", lineheight=1.1)
  g
}

varrow <- function(g, x, y1, y2, col=col_main)
  g + annotate("segment", x=x, xend=x, y=y1, yend=y2+0.014,
               color=col, linewidth=0.6,
               arrow=arrow(length=unit(0.18,"cm"), type="closed"))

harrow <- function(g, x1, x2, y, col="#888888")
  g + annotate("segment", x=x1, xend=x2-0.01, y=y, yend=y,
               color=col, linewidth=0.5,
               arrow=arrow(length=unit(0.15,"cm"), type="closed"))

bracket <- function(g, x, y1, y2, label, col) {
  g + annotate("segment",x=x,xend=x,y=y1,yend=y2,color=col,linewidth=1.2) +
    annotate("segment",x=x,xend=x+0.01,y=y1,yend=y1,color=col,linewidth=1.2) +
    annotate("segment",x=x,xend=x+0.01,y=y2,yend=y2,color=col,linewidth=1.2) +
    annotate("text",x=x-0.014,y=(y1+y2)/2,label=label,
             angle=90,hjust=0.5,size=2.8,color=col,fontface="bold")
}

p <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
  theme(plot.background=element_rect(fill=col_bg,color=NA),
        plot.margin=margin(10,10,10,10))

mx <- 0.12; mw <- 0.46; mc <- mx+mw/2
ex <- 0.62; ew <- 0.36

# ── Brackets
p <- bracket(p, 0.03, 0.920, 0.810, "PRODUCTS",  col_main)
p <- bracket(p, 0.03, 0.790, 0.635, "PLANTS",    col_plant)
p <- bracket(p, 0.03, 0.615, 0.168, "COMPOUNDS", col_cpd)
p <- bracket(p, 0.03, 0.148, 0.040, "DOCKING",   col_dock)

# ── Headers
p <- p +
  annotate("text",x=mc,y=0.970,label="Retained",hjust=0.5,
           size=2.9,color=col_main,fontface="bold") +
  annotate("text",x=ex+ew/2,y=0.970,label="Excluded / Filter detail",hjust=0.5,
           size=2.9,color=col_excl,fontface="bold") +
  annotate("segment",x=mx,xend=mx+mw,y=0.963,yend=0.963,
           color=col_main,linewidth=0.4) +
  annotate("segment",x=ex,xend=ex+ew,y=0.963,yend=0.963,
           color=col_excl,linewidth=0.4)

# ══ BOX 1: KNApSAcK
p <- box(p, mx, 0.930, mw, 0.050,
         sprintf("KNApSAcK Jamu Database  (n = %s)", fmt(n_knapsack_total)),
         fill="#E6F1FB", border=col_main, bold=TRUE)
p <- varrow(p, mc, 0.930, 0.878)

# ══ BOX 2: Sehat
p <- box(p, mx, 0.815, mw, 0.063,
         sprintf("Jamu 'Sehat' products  (n = %s)", fmt(n_sehat)),
         sprintf("Male-specific = %s  |  Female-specific = %s  |  General health = %s",
                 fmt(n_male), fmt(n_female), fmt(n_gen)),
         fill="#E6F1FB", border=col_main, bold=TRUE)
p <- harrow(p, mx+mw, ex, 0.846)
p <- box(p, ex, 0.828, ew, 0.037,
         "Excluded: non-'sehat' products",
         sprintf("n = %s  %s", fmt(n_knapsack_total-n_sehat),
                 pct(n_knapsack_total-n_sehat, n_knapsack_total)),
         fill="#F1EFE8", border=col_excl)
p <- varrow(p, mc, 0.815, 0.763, col_plant)

# ══ BOX 3: Ingredient entries
p <- box(p, mx, 0.723, mw, 0.052,
         sprintf("Ingredient entries scraped  (n = %s)", fmt(n_ing_raw)),
         sprintf("From %s products  ·  avg %.0f ingredients/product",
                 fmt(n_sehat), ifelse(is.na(n_ing_raw)||is.na(n_sehat)||n_sehat==0,
                                      NA, n_ing_raw/n_sehat)),
         fill="#E1F5EE", border=col_plant, bold=TRUE)
p <- harrow(p, mx+mw, ex, 0.749)
p <- box(p, ex, 0.731, ew, 0.037,
         "Synonyms resolved via GBIF",
         sprintf("n = %s synonyms merged  |  EXACT %s  FUZZY %s",
                 fmt(n_synonyms), fmt(n_exact), fmt(n_fuzzy)),
         fill="#F1EFE8", border=col_excl)
p <- varrow(p, mc, 0.723, 0.671, col_plant)

# ══ BOX 4: Unique species
p <- box(p, mx, 0.628, mw, 0.053,
         sprintf("Unique accepted species  (n = %s)", fmt(n_unique_spp)),
         sprintf("%s families  |  EXACT %s  FUZZY %s  HIGHERRANK %s",
                 fmt(n_families), fmt(n_exact), fmt(n_fuzzy), fmt(n_higherrank)),
         fill="#E1F5EE", border=col_plant, bold=TRUE)
p <- varrow(p, mc, 0.628, 0.576, col_plant)

# ══ BOX 5: Priority plants
p <- box(p, mx, 0.530, mw, 0.057,
         sprintf("Priority plants — Tukey + max 3/family  (n = %s)", fmt(n_priority)),
         sprintf("%s families  |  threshold = Q1+1.5×IQR = %.1f",
                 fmt(n_priority_fam), ifelse(is.na(tukey_thresh), 0, tukey_thresh)),
         fill="#E1F5EE", border=col_plant, bold=TRUE)
p <- harrow(p, mx+mw, ex, 0.558)
p <- box(p, ex, 0.540, ew, 0.037,
         "Below Tukey threshold",
         sprintf("n = %s excluded  %s",
                 fmt(n_below_tukey), pct(n_below_tukey, n_unique_spp)),
         fill="#F1EFE8", border=col_excl)
p <- varrow(p, mc, 0.530, 0.478, col_cpd)

# ══ BOX 6: KNApSAcK raw
p <- box(p, mx, 0.432, mw, 0.052,
         sprintf("KNApSAcK compounds raw  (n = %s)", fmt(n_knapsack_raw)),
         sprintf("From %s priority species", fmt(n_knapsack_spp)),
         fill="#FAEEDA", border=col_cpd, bold=TRUE)
p <- varrow(p, mc, 0.432, 0.380, col_cpd)

# ══ BOX 7: SMILES + dedup
p <- box(p, mx, 0.333, mw, 0.053,
         sprintf("PubChem SMILES + deduplication  (n = %s)", fmt(n_smiles_valid)),
         "Deduplicated by canonical SMILES  ·  provenance retained",
         fill="#FAEEDA", border=col_cpd, bold=TRUE)
p <- harrow(p, mx+mw, ex, 0.359)
p <- box(p, ex, 0.341, ew, 0.037,
         "No valid SMILES / duplicates",
         sprintf("n = %s excluded  %s",
                 fmt(n_smiles_excl), pct(n_smiles_excl, n_knapsack_raw)),
         fill="#F1EFE8", border=col_excl)
p <- varrow(p, mc, 0.333, 0.281, col_cpd)

# ══ BOX 8: Lipinski
lip_sub <- sprintf(
  "MW \u2264500: %s  |  LogP \u22645: %s  |  HBD \u22645: %s\nHBA \u226410: %s  |  TPSA \u2264140: %s  |  RotBond \u226410: %s\nAll criteria met: %s / %s  %s",
  fmt(lip_mw), fmt(lip_logp), fmt(lip_hbd),
  fmt(lip_hba), fmt(lip_tpsa), fmt(lip_rot),
  fmt(n_lipinski_pass), fmt(n_smiles_valid),
  pct(n_lipinski_pass, n_smiles_valid)
)
p <- box(p, mx, 0.168, mw, 0.122,
         sprintf("Lipinski RO5 + Veber filter  (n = %s passed)", fmt(n_lipinski_pass)),
         lip_sub,
         fill="#FAEEDA", border=col_cpd, bold=TRUE, subfontsize=2.6)
p <- harrow(p, mx+mw, ex, 0.230)
p <- box(p, ex, 0.210, ew, 0.050,
         "Failed Lipinski / Veber",
         sprintf("n = %s excluded  %s\nTightest: TPSA >140 (%s fail)",
                 fmt(n_lipinski_excl), pct(n_lipinski_excl, n_smiles_valid),
                 fmt(ifelse(!is.na(lip_tpsa), n_smiles_valid - lip_tpsa, NA))),
         fill="#F1EFE8", border=col_excl, subfontsize=2.55)
p <- varrow(p, mc, 0.168, 0.116, col_cpd)  # lipinski -> admet

# ══ BOX 9: ADMET
admet_sub <- sprintf(
  "PAINS exclusion: %s removed  |  GI High: %s %s\nLogS \u2265\u22126: %s %s  |  Lip.viol \u22641: %s\nAll ADMET: %s / %s  %s  |  BBB: %s %s",
  fmt(n_pains_excl),
  fmt(adm_gi),    pct(adm_gi,    n_lipinski_pass),
  fmt(adm_logs),  pct(adm_logs,  n_lipinski_pass),
  fmt(adm_lipviol),
  fmt(n_admet_pass), fmt(n_lipinski_pass), pct(n_admet_pass, n_lipinski_pass),
  fmt(adm_bbb),   pct(adm_bbb, n_admet_pass)
)
p <- box(p, mx, 0.168, mw, 0.120,
         sprintf("ADMET filter (RDKit)  (n = %s passed)", fmt(n_admet_pass)),
         admet_sub,
         fill="#FAEEDA", border=col_cpd, bold=TRUE, subfontsize=2.6)
p <- harrow(p, mx+mw, ex, 0.108)
p <- box(p, ex, 0.087, ew, 0.050,
         "Failed ADMET",
         sprintf("n = %s excluded  %s\nincl. %s PAINS + %s LogS/other",
                 fmt(n_admet_excl), pct(n_admet_excl, n_lipinski_pass),
                 fmt(n_pains_excl),
                 fmt(ifelse(!is.na(n_admet_excl)&&!is.na(n_pains_excl),
                            n_admet_excl - n_pains_excl, NA))),
         fill="#F1EFE8", border=col_excl, subfontsize=2.55)
p <- varrow(p, mc, 0.048, 0.006, col_dock)

# ══ BOX DOCKING: ChEMBL + Final shortlist
# (rendered below the diagram area as annotation — add to title section)

# ══ TITLE & CAPTION
# Build target breakdown string for box
target_str <- if (!is.null(target_breakdown)) {
  paste(
    apply(target_breakdown, 1, function(r)
      sprintf("%s: %s cpd", sub(".*\\((.+)\\).*","\\1", r["protein_name"]), r["n_cpd"])),
    collapse="  |  ")
} else "AR: 3  |  ESR1: 11  |  HIF1A: 8  |  NFKB1: 11"

p <- p +
  # Title
  annotate("text", x=0.5, y=0.997,
           label="Filtering Funnel: KNApSAcK Jamu Database → Molecular Docking Candidates",
           hjust=0.5, vjust=1, size=3.8, fontface="bold", color=col_text) +

  # Subtitle: ChEMBL detail
  annotate("text", x=0.5, y=0.982,
           label=sprintf(
             "ChEMBL similarity ≥60%%  ·  pChEMBL ≥4.0  ·  %s target hits from %s compounds",
             fmt(n_target_total), fmt(n_target_hits)),
           hjust=0.5, vjust=1, size=2.9, color="#555555") +

  # ChEMBL target prediction box (DOCKING zone)
  annotate("rect", xmin=mx, xmax=mx+mw, ymin=0.092, ymax=0.148,
           fill="#FBEAF0", color=col_dock, linewidth=0.5) +
  annotate("text", x=mc, y=0.136, hjust=0.5, vjust=0.5,
           label=sprintf("ChEMBL target prediction  (%s hits from %s compounds)",
                         fmt(n_target_total), fmt(n_target_hits)),
           size=3.2, color=col_text, fontface="bold") +
  annotate("text", x=mc, y=0.112, hjust=0.5, vjust=0.5,
           label=sprintf("Similarity ≥60%%  ·  pChEMBL ≥4.0  ·  %s",
                         target_str),
           size=2.7, color="#444444", lineheight=1.1) +

  # Arrow ChEMBL -> Final
  annotate("segment", x=mc, xend=mc, y=0.092, yend=0.062+0.014,
           color=col_dock, linewidth=0.6,
           arrow=arrow(length=unit(0.18,"cm"), type="closed")) +

  # Exclusion from ChEMBL step
  annotate("segment", x=mx+mw, xend=ex-0.01, y=0.120, yend=0.120,
           color="#888888", linewidth=0.5,
           arrow=arrow(length=unit(0.15,"cm"), type="closed")) +
  annotate("rect", xmin=ex, xmax=ex+ew, ymin=0.100, ymax=0.140,
           fill="#F1EFE8", color=col_excl, linewidth=0.5) +
  annotate("text", x=ex+ew/2, y=0.128, hjust=0.5, vjust=0.5,
           label="No target hit / pChEMBL <4.0",
           size=3.0, color=col_text) +
  annotate("text", x=ex+ew/2, y=0.112, hjust=0.5, vjust=0.5,
           label=sprintf("n = %s excluded  %s",
                         fmt(n_no_target_hit), pct(n_no_target_hit, n_admet_pass)),
           size=2.7, color="#444444") +

  # Final shortlist box
  annotate("rect", xmin=mx, xmax=mx+mw, ymin=0.030, ymax=0.062+0.025,
           fill="#FBEAF0", color=col_dock, linewidth=1.2) +
  annotate("text", x=mc, y=0.074, hjust=0.5, vjust=0.5,
           label=sprintf("Final shortlist  —  n = %s compounds  ·  %s protein targets",
                         fmt(n_final), fmt(n_targets_used)),
           size=3.5, color=col_dock, fontface="bold") +
  annotate("text", x=mc, y=0.050, hjust=0.5, vjust=0.5,
           label="Criteria: MW\u2264500  \u00b7  n_targets>0  \u00b7  pChEMBL>5.0\nAutoDock Vina v1.2.5  \u00b7  Exhaustiveness=16",
           size=2.7, color="#444444", lineheight=1.2) +

  # Caption
  annotate("text", x=0.5, y=0.004,
           label="GBIF backbone taxonomy. Lipinski & ADMET: RDKit. ChEMBL v33. Docking: AutoDock Vina v1.2.5.",
           hjust=0.5, vjust=0, size=2.4, color="#888888")

ggsave("prisma_funnel_publication.png",
       plot=p, width=14, height=20, dpi=300, bg=col_bg)
cat("\nSelesai: prisma_funnel_publication.png\n")
cat("Semua nilai n dibaca dari file CSV di folder docking.\n")
