# ============================================================
# Jamu Formulation Score (JFS)
# Menggabungkan ethnopharmacology tendency + molecular target
# untuk merekomendasikan formulasi jamu Male / Female / General
# Target: AR (P10275), ESR1 (P03372), HIF1A (Q16665), NFKB1 (P19838)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(boot)

set.seed(42)

# ════════════════════════════════════════════════════════════
# 1. BACA DATA
# ════════════════════════════════════════════════════════════

plants  <- read.csv("bobot_tanaman.csv",      stringsAsFactors = FALSE)
targets <- read.csv("target_predictions.csv", stringsAsFactors = FALSE)
cpd_raw <- read.csv("compounds_final.csv",    stringsAsFactors = FALSE)

# Rename jika perlu
if ("SPECIES" %in% names(plants))
  names(plants)[names(plants) == "SPECIES"] <- "TANAMAN"

cat("Plants loaded    :", nrow(plants), "\n")
cat("Target pred      :", nrow(targets), "\n")

# ════════════════════════════════════════════════════════════
# 2. FILTER: hanya 4 target yang digunakan, exclude Se-Methyl
# ════════════════════════════════════════════════════════════

TARGET_IDS <- c("P10275", "P03372", "Q16665", "P19838")
# P10275 = AR, P03372 = ESR1, Q16665 = HIF1A, P19838 = NFKB1

tgt <- targets %>%
  filter(
    uniprot_id %in% TARGET_IDS,
    compound_name != "Se-Methylselenocysteine"
  ) %>%
  # Deduplikasi: ambil pChEMBL tertinggi per compound-protein pair
  group_by(compound_name, species, uniprot_id, protein_name) %>%
  summarise(pchembl = max(pchembl_value, na.rm = TRUE),
            .groups = "drop")

cat("\nTarget hits (4 proteins, dedup):", nrow(tgt), "\n")
print(tgt %>% select(compound_name, species, protein_name, pchembl) %>%
        arrange(protein_name, desc(pchembl)))

# ════════════════════════════════════════════════════════════
# 3. TARGET RELEVANCE MATRIX
# Nilai 0–1 berdasarkan relevansi biologis terhadap kategori jamu
# Referensi: fisiologi olahraga & endokrinologi
# ════════════════════════════════════════════════════════════
#
#  AR    → Male hormonal axis (testosterone/DHT receptor)
#           → Male stamina tinggi, Female rendah, General sedang
#  ESR1  → Female hormonal axis (estrogen receptor α)
#           → Female tinggi, Male rendah, General sedang
#  HIF1A → Hypoxia-inducible factor: endurance, VO2max
#           → Universal untuk stamina semua kategori
#  NFKB1 → NF-κB: anti-inflamasi & recovery pasca latihan
#           → Universal, sedikit lebih relevan Male (DOMS recovery)

trm <- data.frame(
  uniprot_id   = c("P10275","P03372","Q16665","P19838"),
  protein_short= c("AR",    "ESR1",  "HIF1A", "NFKB1"),
  rel_male     = c(1.00,    0.10,    0.70,    0.80),
  rel_female   = c(0.10,    1.00,    0.70,    0.70),
  rel_general  = c(0.40,    0.40,    1.00,    0.90),
  justification= c(
    "Androgen receptor: testosteron/DHT, massa otot, stamina pria",
    "Estrogen receptor alpha: hormonal wanita, bone density, endurance",
    "HIF-1alpha: adaptasi hipoksia, VO2max, ketahanan aerobik",
    "NF-kB1: anti-inflamasi, DOMS recovery, imunitas latihan"
  ),
  stringsAsFactors = FALSE
)

cat("\n=== Target Relevance Matrix ===\n")
print(trm %>% select(protein_short, rel_male, rel_female, rel_general,
                     justification))

# ════════════════════════════════════════════════════════════
# 4. HITUNG BIOACTIVITY SCORE PER TANAMAN PER KATEGORI
# ════════════════════════════════════════════════════════════
# BAS(plant, cat) = Σ [ norm_pChEMBL(cpd) × rel(target, cat) ]
#                   untuk semua cpd dari tanaman tsb
#
# norm_pChEMBL: skala 0–1, anchor di IC50 = 10µM (pChEMBL=5)
#   norm = (pchembl - 4) / (max_pchembl - 4)
#   di mana 4 = batas minimal ChEMBL, max = max observed

max_pchembl <- max(tgt$pchembl)
cat(sprintf("\npChEMBL range: %.2f – %.2f\n", min(tgt$pchembl), max_pchembl))

tgt <- tgt %>%
  mutate(norm_pchembl = (pchembl - 4) / (max_pchembl - 4))

# Gabung dengan target relevance
tgt_full <- tgt %>%
  left_join(trm %>% select(uniprot_id, protein_short,
                            rel_male, rel_female, rel_general),
            by = "uniprot_id") %>%
  mutate(
    contrib_male    = norm_pchembl * rel_male,
    contrib_female  = norm_pchembl * rel_female,
    contrib_general = norm_pchembl * rel_general
  )

# Agregasi per tanaman
bas <- tgt_full %>%
  group_by(species) %>%
  summarise(
    BAS_male    = sum(contrib_male,    na.rm = TRUE),
    BAS_female  = sum(contrib_female,  na.rm = TRUE),
    BAS_general = sum(contrib_general, na.rm = TRUE),
    n_compounds = n_distinct(compound_name),
    n_targets   = n_distinct(uniprot_id),
    compounds   = paste(unique(compound_name), collapse = "; "),
    proteins    = paste(unique(protein_short), collapse = "/"),
    .groups     = "drop"
  ) %>%
  rename(TANAMAN = species)

cat("\n=== Bioactivity Score per tanaman ===\n")
print(bas %>% select(TANAMAN, BAS_male, BAS_female, BAS_general,
                     n_compounds, proteins) %>%
        arrange(desc(BAS_general)))

# ════════════════════════════════════════════════════════════
# 5. HITUNG TENDENCY SCORE PER TANAMAN PER KATEGORI
# ════════════════════════════════════════════════════════════
# TS(plant, cat) = FREQ_cat / total_freq   (proporsi relatif)
# Tanaman dengan total_freq = 0 → TS = 0 (tidak digunakan di jamu)

plants <- plants %>%
  mutate(
    total_freq = MALE_FREQ + FEM_FREQ + NEU_FREQ,
    TS_male    = ifelse(total_freq > 0, MALE_FREQ / total_freq, 0),
    TS_female  = ifelse(total_freq > 0, FEM_FREQ  / total_freq, 0),
    TS_general = ifelse(total_freq > 0, NEU_FREQ  / total_freq, 0),
    # Mixed tanaman yang tidak spesifik ke satu kategori
    # mendapat score dari semua frekuensi yg dimiliki
    TS_mixed   = ifelse(total_freq > 0, total_freq / max(total_freq), 0)
  )

cat("\n=== Tendency Score (5 tanaman teratas per kategori) ===\n")
cat("Male:\n")
plants %>% arrange(desc(TS_male)) %>% head(5) %>%
  select(TANAMAN, MALE_FREQ, FEM_FREQ, NEU_FREQ, TS_male) %>% print()
cat("Female:\n")
plants %>% arrange(desc(TS_female)) %>% head(5) %>%
  select(TANAMAN, MALE_FREQ, FEM_FREQ, NEU_FREQ, TS_female) %>% print()

# ════════════════════════════════════════════════════════════
# 6. HITUNG JAMU FORMULATION SCORE (JFS)
# ════════════════════════════════════════════════════════════
# JFS(plant, cat) = TS(plant, cat) × BAS(plant, cat) × final_score_norm
#
# final_score_norm = final_score / max(final_score)
# → memperhitungkan seberapa sering tanaman muncul di jamu
#
# Untuk tanaman TANPA data senyawa (BAS = 0):
#   JFS_ethnopharm = final_score_norm × TS   (fallback, BAS = 1)
#   Ini berarti tanaman tersebut dinilai dari ethnopharmacology saja

# Normalize final_score
plants <- plants %>%
  mutate(fs_norm = final_score / max(final_score))

# Merge BAS ke plants
jfs <- plants %>%
  left_join(bas, by = "TANAMAN") %>%
  # Tanaman tanpa BAS (no compound data) → BAS = 0 tapi tetap masuk
  mutate(
    BAS_male    = ifelse(is.na(BAS_male),    0, BAS_male),
    BAS_female  = ifelse(is.na(BAS_female),  0, BAS_female),
    BAS_general = ifelse(is.na(BAS_general), 0, BAS_general),
    n_compounds = ifelse(is.na(n_compounds), 0, n_compounds),
    has_compound = n_compounds > 0,

    # JFS: gabungan tendency + bioactivity + ethnopharm weight
    # Untuk tanaman dengan senyawa:
    #   JFS = TS × BAS × fs_norm (semua faktor berkontribusi)
    # Untuk tanaman tanpa senyawa:
    #   JFS = TS × fs_norm × 0.3 (penalty krn tidak ada data bioaktivitas)
    JFS_male = ifelse(
      has_compound,
      TS_male    * BAS_male    * fs_norm,
      TS_male    * fs_norm     * 0.3
    ),
    JFS_female = ifelse(
      has_compound,
      TS_female  * BAS_female  * fs_norm,
      TS_female  * fs_norm     * 0.3
    ),
    JFS_general = ifelse(
      has_compound,
      TS_general * BAS_general * fs_norm,
      TS_general * fs_norm     * 0.3
    ),

    # Overall JFS = max dari ketiga kategori (potensi terbaik tanaman)
    JFS_overall = pmax(JFS_male, JFS_female, JFS_general)
  ) %>%
  select(TANAMAN, FAMILI, group, total_freq, final_score, fs_norm,
         TS_male, TS_female, TS_general,
         BAS_male, BAS_female, BAS_general,
         JFS_male, JFS_female, JFS_general, JFS_overall,
         has_compound, n_compounds, n_targets, compounds, proteins)

cat("\n=== Jamu Formulation Score — semua tanaman ===\n")
print(jfs %>% filter(total_freq > 0) %>%
        select(TANAMAN, group, JFS_male, JFS_female, JFS_general, has_compound) %>%
        arrange(desc(JFS_overall <- pmax(JFS_male, JFS_female, JFS_general))))

# ════════════════════════════════════════════════════════════
# 7. UJI STATISTIK
# ════════════════════════════════════════════════════════════

cat("\n\n═══════════════════════════════════════════════════════\n")
cat("STATISTICAL TESTING\n")
cat("═══════════════════════════════════════════════════════\n")

active_plants <- jfs %>% filter(total_freq > 0)

# ── 7a. Wilcoxon rank-sum: apakah distribusi JFS berbeda antar kategori ──
cat("\n[1] Wilcoxon rank-sum test: Male vs Female JFS\n")
w_mf <- wilcox.test(active_plants$JFS_male, active_plants$JFS_female,
                    paired = FALSE, exact = FALSE)
cat(sprintf("    W = %.1f, p = %.4f %s\n", w_mf$statistic, w_mf$p.value,
            ifelse(w_mf$p.value < 0.05, "(significant)", "(not significant)")))

cat("\n[2] Wilcoxon rank-sum test: Male vs General JFS\n")
w_mg <- wilcox.test(active_plants$JFS_male, active_plants$JFS_general,
                    paired = FALSE, exact = FALSE)
cat(sprintf("    W = %.1f, p = %.4f %s\n", w_mg$statistic, w_mg$p.value,
            ifelse(w_mg$p.value < 0.05, "(significant)", "(not significant)")))

cat("\n[3] Wilcoxon rank-sum test: Female vs General JFS\n")
w_fg <- wilcox.test(active_plants$JFS_female, active_plants$JFS_general,
                    paired = FALSE, exact = FALSE)
cat(sprintf("    W = %.1f, p = %.4f %s\n", w_fg$statistic, w_fg$p.value,
            ifelse(w_fg$p.value < 0.05, "(significant)", "(not significant)")))

# ── 7b. Kruskal-Wallis: overall difference across 3 categories ──
cat("\n[4] Kruskal-Wallis test: perbedaan JFS antar 3 kategori\n")
jfs_long <- active_plants %>%
  select(TANAMAN, JFS_male, JFS_female, JFS_general) %>%
  pivot_longer(cols = starts_with("JFS"),
               names_to = "category", values_to = "JFS") %>%
  mutate(category = sub("JFS_", "", category))

kw <- kruskal.test(JFS ~ category, data = jfs_long)
cat(sprintf("    chi-sq = %.3f, df = %d, p = %.4f %s\n",
            kw$statistic, kw$parameter, kw$p.value,
            ifelse(kw$p.value < 0.05, "(significant)", "(not significant)")))

# ── 7c. Permutation test: top-5 tanaman per kategori ──
cat("\n[5] Permutation test: apakah top-5 tanaman signifikan?\n")
N_PERM <- 10000

for (cat_name in c("male", "female", "general")) {
  jfs_col   <- paste0("JFS_", cat_name)
  top5_obs  <- active_plants %>%
    arrange(desc(.data[[jfs_col]])) %>%
    slice_head(n = 5) %>%
    pull(.data[[jfs_col]])
  obs_score <- sum(top5_obs)

  null_dist <- replicate(N_PERM, {
    sampled <- sample(active_plants[[jfs_col]], size = 5, replace = FALSE)
    sum(sampled)
  })

  p_perm <- mean(null_dist >= obs_score)
  ci_lo  <- quantile(null_dist, 0.025)
  ci_hi  <- quantile(null_dist, 0.975)

  cat(sprintf("    %s: observed=%.4f | null CI=[%.4f, %.4f] | p=%.4f %s\n",
              toupper(cat_name), obs_score, ci_lo, ci_hi, p_perm,
              ifelse(p_perm < 0.05, "✓ SIGNIFICANT", "✗ not significant")))
}

# ── 7d. Bootstrap CI untuk JFS top tanaman ──
cat("\n[6] Bootstrap 95% CI untuk JFS tanaman teratas (B=2000)\n")

for (cat_name in c("male", "female", "general")) {
  jfs_col <- paste0("JFS_", cat_name)
  top_plant <- active_plants %>%
    arrange(desc(.data[[jfs_col]])) %>%
    slice_head(n = 1)

  # Bootstrap CI via resampling dari pool tanaman yang punya freq > 0
  vals <- active_plants[[jfs_col]]
  boot_fn <- function(data, i) mean(data[i])
  boot_res <- boot(data = vals, statistic = boot_fn, R = 2000)
  ci <- boot.ci(boot_res, type = "perc", conf = 0.95)

  cat(sprintf("    %s top: %s | JFS=%.4f | 95%% CI=[%.4f, %.4f]\n",
              toupper(cat_name),
              top_plant$TANAMAN,
              top_plant[[jfs_col]],
              ci$percent[4], ci$percent[5]))
}

# ════════════════════════════════════════════════════════════
# 8. RANKING & REKOMENDASI FORMULASI
# ════════════════════════════════════════════════════════════
cat("\n\n═══════════════════════════════════════════════════════\n")
cat("REKOMENDASI FORMULASI JAMU\n")
cat("═══════════════════════════════════════════════════════\n")

# Fungsi seleksi: top-N dengan max 2 per famili
select_formulation <- function(data, jfs_col, n = 5, max_per_family = 2) {
  data %>%
    filter(total_freq > 0) %>%
    arrange(desc(.data[[jfs_col]])) %>%
    group_by(FAMILI) %>%
    mutate(rank_in_family = row_number()) %>%
    ungroup() %>%
    filter(rank_in_family <= max_per_family) %>%
    slice_head(n = n) %>%
    select(TANAMAN, FAMILI, group,
           TS    = !!sym(sub("JFS","TS", jfs_col)),
           BAS   = !!sym(sub("JFS","BAS", jfs_col)),
           JFS   = !!sym(jfs_col),
           has_compound, compounds, proteins) %>%
    mutate(rank = row_number())
}

form_male    <- select_formulation(jfs, "JFS_male",    n = 5)
form_female  <- select_formulation(jfs, "JFS_female",  n = 5)
form_general <- select_formulation(jfs, "JFS_general", n = 5)

cat("\n── FORMULASI MALE-SPECIFIC (stamina pria) ──────────────\n")
cat("Target dominan: AR (testosterone/DHT) + NFKB1 (recovery)\n\n")
print(form_male %>% select(rank, TANAMAN, FAMILI, TS, BAS, JFS,
                            has_compound, proteins))

cat("\n── FORMULASI FEMALE-SPECIFIC (stamina wanita) ──────────\n")
cat("Target dominan: ESR1 (estrogen) + HIF1A (endurance)\n\n")
print(form_female %>% select(rank, TANAMAN, FAMILI, TS, BAS, JFS,
                              has_compound, proteins))

cat("\n── FORMULASI GENERAL HEALTH (stamina universal) ────────\n")
cat("Target dominan: HIF1A (VO2max) + NFKB1 (anti-inflamasi)\n\n")
print(form_general %>% select(rank, TANAMAN, FAMILI, TS, BAS, JFS,
                               has_compound, proteins))

# ════════════════════════════════════════════════════════════
# 9. SIMPAN OUTPUT
# ════════════════════════════════════════════════════════════
write.csv(jfs, "jamu_formulation_scores.csv", row.names = FALSE)

form_all <- bind_rows(
  form_male    %>% mutate(formulation = "Male-specific"),
  form_female  %>% mutate(formulation = "Female-specific"),
  form_general %>% mutate(formulation = "General health")
)
write.csv(form_all, "jamu_formulation_recommended.csv", row.names = FALSE)

write.csv(trm %>% select(protein_short, rel_male, rel_female,
                          rel_general, justification),
          "target_relevance_matrix.csv", row.names = FALSE)

cat("\nFile disimpan:\n")
cat("  jamu_formulation_scores.csv      — JFS semua 70 tanaman\n")
cat("  jamu_formulation_recommended.csv — top-5 per kategori\n")
cat("  target_relevance_matrix.csv      — target relevance matrix\n")

# ════════════════════════════════════════════════════════════
# 10. VISUALISASI
# ════════════════════════════════════════════════════════════

# ── Plot 1: Heatmap JFS semua tanaman aktif ──────────────────
jfs_plot <- jfs %>%
  filter(total_freq > 0) %>%
  select(TANAMAN, JFS_male, JFS_female, JFS_general) %>%
  pivot_longer(cols = starts_with("JFS"),
               names_to = "category",
               values_to = "JFS") %>%
  mutate(
    category = recode(category,
                      JFS_male    = "Male-specific",
                      JFS_female  = "Female-specific",
                      JFS_general = "General health"),
    TANAMAN_short = sub("(\\S+) (\\S).*", "\\1 \\2.", TANAMAN)
  )

# Order tanaman by overall JFS
plant_order <- jfs %>%
  filter(total_freq > 0) %>%
  mutate(JFS_max = pmax(JFS_male, JFS_female, JFS_general)) %>%
  arrange(JFS_max) %>%
  pull(TANAMAN) %>%
  sub("(\\S+) (\\S).*", "\\1 \\2.", .)

jfs_plot <- jfs_plot %>%
  mutate(TANAMAN_short = factor(TANAMAN_short, levels = plant_order))

p_heat <- ggplot(jfs_plot, aes(x = category, y = TANAMAN_short, fill = JFS)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(JFS > 0.001, sprintf("%.3f", JFS), "")),
            size = 2.5, color = "white") +
  scale_fill_gradientn(
    colors  = c("#F1EFE8","#FAC775","#BA7517","#633806","#412402"),
    name    = "JFS",
    limits  = c(0, NA)
  ) +
  scale_x_discrete(position = "top") +
  labs(
    title    = "Jamu Formulation Score (JFS) per tanaman",
    subtitle = "Kombinasi ethnopharmacology tendency × bioactivity (AR, ESR1, HIF1A, NFKB1)",
    x = NULL, y = NULL,
    caption  = "JFS = TS × BAS × fs_norm  |  Tanaman tanpa data senyawa: JFS × 0.3 penalty"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, color = "grey40"),
    axis.text.y     = element_text(face = "italic", size = 8.5),
    axis.text.x     = element_text(face = "bold", size = 10),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid      = element_blank()
  )

# ── Plot 2: Bar chart top-5 per kategori ─────────────────────
form_plot <- form_all %>%
  mutate(
    TANAMAN_short = sub("(\\S+) (\\S).*", "\\1 \\2.", TANAMAN),
    formulation   = factor(formulation,
                           levels = c("Male-specific","Female-specific","General health")),
    bar_label     = sprintf("%.3f", JFS),
    has_cpd_label = ifelse(has_compound, "*", "")
  )

col_form <- c(
  "Male-specific"  = "#1565A0",
  "Female-specific"= "#993556",
  "General health" = "#0F6E56"
)

p_bar <- ggplot(form_plot,
                aes(x = reorder(TANAMAN_short, JFS), y = JFS,
                    fill = formulation)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(bar_label, has_cpd_label)),
            hjust = -0.1, size = 3.0, color = "grey20") +
  coord_flip() +
  facet_wrap(~ formulation, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = col_form, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(
    title    = "Top-5 tanaman per formulasi jamu",
    subtitle = "* = tanaman dengan data senyawa tervalidasi  |  max 2 tanaman per famili",
    x = NULL, y = "Jamu Formulation Score (JFS)",
    caption  = "Seleksi: Tukey threshold + max 2 per famili"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 8.5, color = "grey40"),
    strip.text      = element_text(face = "bold", size = 10),
    axis.text.y     = element_text(face = "italic", size = 9),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major.y = element_blank()
  )

# ── Plot 3: Scatter TS vs BAS per kategori ───────────────────
jfs_scatter <- jfs %>%
  filter(total_freq > 0) %>%
  pivot_longer(
    cols      = c(JFS_male, JFS_female, JFS_general),
    names_to  = "cat", values_to = "JFS"
  ) %>%
  mutate(
    TS  = case_when(cat=="JFS_male"    ~ TS_male,
                    cat=="JFS_female"  ~ TS_female,
                    cat=="JFS_general" ~ TS_general),
    BAS = case_when(cat=="JFS_male"    ~ BAS_male,
                    cat=="JFS_female"  ~ BAS_female,
                    cat=="JFS_general" ~ BAS_general),
    category = recode(cat,
                      JFS_male    = "Male-specific",
                      JFS_female  = "Female-specific",
                      JFS_general = "General health"),
    TANAMAN_short = sub("(\\S+) (\\S).*", "\\1 \\2.", TANAMAN),
    top5 = TANAMAN %in% form_all$TANAMAN
  )

p_scatter <- ggplot(jfs_scatter,
                    aes(x = TS, y = BAS, size = JFS,
                        color = category, alpha = top5)) +
  geom_point(stroke = 0.5, shape = 21,
             aes(fill = category), color = "white") +
  geom_text(data = jfs_scatter %>% filter(top5 & JFS > 0),
            aes(label = TANAMAN_short),
            size = 2.4, nudge_y = 0.02, fontface = "italic",
            show.legend = FALSE) +
  scale_fill_manual(values  = col_form, name = "Formulasi") +
  scale_color_manual(values = col_form, guide = "none") +
  scale_alpha_manual(values = c("TRUE"=1.0,"FALSE"=0.3), guide="none") +
  scale_size_continuous(range = c(1.5, 8), name = "JFS") +
  facet_wrap(~ category, ncol = 3) +
  labs(
    title    = "Tendency Score vs Bioactivity Score",
    subtitle = "Ukuran titik ∝ JFS  |  Label = tanaman terpilih dalam formulasi",
    x = "Tendency Score (TS) — proporsi frekuensi jamu",
    y = "Bioactivity Score (BAS) — sum pChEMBL × target relevance"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 8.5, color = "grey40"),
    strip.text      = element_text(face = "bold", size = 10),
    plot.background = element_rect(fill = "white", color = NA)
  )

# ── Simpan plots ─────────────────────────────────────────────
ggsave("jfs_heatmap.png",  plot = p_heat,    width=10, height=12, dpi=300, bg="white")
ggsave("jfs_top5_bar.png", plot = p_bar,     width=10, height=13, dpi=300, bg="white")
ggsave("jfs_scatter.png",  plot = p_scatter, width=14, height=6,  dpi=300, bg="white")

cat("\nPlot disimpan:\n")
cat("  jfs_heatmap.png   — heatmap JFS semua tanaman aktif\n")
cat("  jfs_top5_bar.png  — bar chart top-5 per formulasi\n")
cat("  jfs_scatter.png   — scatter TS vs BAS\n")
cat("\nSelesai.\n")
