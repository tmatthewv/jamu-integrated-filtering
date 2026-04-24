# ============================================================
# Filtering Funnel Report — Publication Quality
# ggplot2 + patchwork: main funnel + detail panels
# ============================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# ── Palette ──────────────────────────────────────────────────
col_prod  <- "#1565A0"
col_plant <- "#0F6E56"
col_cpd   <- "#BA7517"
col_dock  <- "#993556"
col_excl  <- "#888888"
col_bg    <- "#FFFFFF"
col_text  <- "#1A1A1A"

# ═══════════════════════════════════════════════════════════
# PANEL 1 — Main funnel
# ═══════════════════════════════════════════════════════════

funnel <- data.frame(
  step  = factor(c(
    "KNApSAcK Jamu\n(5,310 products)",
    "Jamu 'Sehat'\n(368 products)",
    "Ingredient entries\n(15,422 entries)",
    "Unique species\npost-GBIF (234 spp.)",
    "Priority plants\nTukey + max 3/family\n(41 plants)",
    "KNApSAcK\ncompounds raw\n(3,699)",
    "PubChem SMILES\n+ dedup (2,246)",
    "Lipinski RO5\n+ Veber (1,464)",
    "ADMET filter\n(1,292)",
    "ChEMBL target\nprediction",
    "Final shortlist\n8 compounds\n6 targets"
  ), levels = rev(c(
    "KNApSAcK Jamu\n(5,310 products)",
    "Jamu 'Sehat'\n(368 products)",
    "Ingredient entries\n(15,422 entries)",
    "Unique species\npost-GBIF (234 spp.)",
    "Priority plants\nTukey + max 3/family\n(41 plants)",
    "KNApSAcK\ncompounds raw\n(3,699)",
    "PubChem SMILES\n+ dedup (2,246)",
    "Lipinski RO5\n+ Veber (1,464)",
    "ADMET filter\n(1,292)",
    "ChEMBL target\nprediction",
    "Final shortlist\n8 compounds\n6 targets"
  ))),
  n     = c(5310, 368, 15422, 234, 41, 3699, 2246, 1464, 1292, NA, 8),
  phase = c("Products","Products","Plants","Plants","Plants",
            "Compounds","Compounds","Compounds","Compounds","Compounds","Docking")
)

phase_colors <- c(
  Products  = col_prod,
  Plants    = col_plant,
  Compounds = col_cpd,
  Docking   = col_dock
)

# Bar width proportional to log scale for readability
funnel$bar_n  <- ifelse(is.na(funnel$n), 8, funnel$n)
funnel$log_n  <- log10(pmax(funnel$bar_n, 1))
funnel$max_log <- max(funnel$log_n, na.rm = TRUE)

p1 <- ggplot(funnel, aes(y = step, x = log_n, fill = phase)) +
  geom_col(width = 0.68, color = "white", linewidth = 0.3) +
  geom_text(
    aes(x = log_n + 0.06,
        label = ifelse(is.na(n), "n = 8", scales::comma(n))),
    hjust = 0, size = 3.2, color = col_text, fontface = "bold"
  ) +
  scale_fill_manual(
    values = phase_colors,
    name   = "Phase",
    guide  = guide_legend(reverse = FALSE)
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.25)),
    breaks = 0:5,
    labels = c("1", "10", "100", "1K", "10K", "100K")
  ) +
  labs(
    x     = "Count (log scale)",
    y     = NULL,
    title = "Filtering funnel: KNApSAcK Jamu → Molecular docking"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(size = 13, face = "bold", color = col_text,
                                    margin = margin(b = 8)),
    axis.text.y      = element_text(size = 9.5, color = col_text, lineheight = 1.2),
    axis.text.x      = element_text(size = 9, color = "#555555"),
    axis.title.x     = element_text(size = 9, color = "#555555"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.4),
    legend.position    = "right",
    legend.title       = element_text(size = 9, face = "bold"),
    legend.text        = element_text(size = 9),
    plot.background    = element_rect(fill = col_bg, color = NA),
    panel.background   = element_rect(fill = col_bg, color = NA)
  )

# ═══════════════════════════════════════════════════════════
# PANEL 2 — Lipinski detail
# ═══════════════════════════════════════════════════════════

lip <- data.frame(
  criterion = factor(c("MW ≤ 500","LogP ≤ 5","HBD ≤ 5","HBA ≤ 10",
                        "TPSA ≤ 140","RotBond ≤ 10","All passed"),
                      levels = rev(c("MW ≤ 500","LogP ≤ 5","HBD ≤ 5","HBA ≤ 10",
                                     "TPSA ≤ 140","RotBond ≤ 10","All passed"))),
  passed   = c(1911, 1948, 1861, 1933, 1817, 1974, 1464),
  excluded = c(335, 298, 385, 313, 429, 272, 782),
  total    = 2246
)
lip$pct_pass <- round(lip$passed / lip$total * 100, 1)

p2 <- ggplot(lip, aes(y = criterion)) +
  geom_col(aes(x = total), fill = "#EEEEEE", width = 0.65) +
  geom_col(aes(x = passed,
               fill = ifelse(criterion == "All passed", col_cpd, "#EF9F27")),
           width = 0.65, show.legend = FALSE) +
  geom_text(aes(x = passed + 30,
                label = paste0(scales::comma(passed), " (", pct_pass, "%)")),
            hjust = 0, size = 3.1, color = col_text, fontface = "bold") +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, 2700),
    breaks = c(0, 500, 1000, 1500, 2000),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    x     = "Compounds passing (out of 2,246)",
    y     = NULL,
    title = "Lipinski RO5 + Veber — per-criterion detail"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title         = element_text(size = 11, face = "bold", color = col_text,
                                      margin = margin(b = 6)),
    axis.text.y        = element_text(size = 9.5, color = col_text),
    axis.text.x        = element_text(size = 8.5, color = "#555555"),
    axis.title.x       = element_text(size = 8.5, color = "#555555"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.4),
    plot.background    = element_rect(fill = col_bg, color = NA),
    panel.background   = element_rect(fill = col_bg, color = NA)
  )

# ═══════════════════════════════════════════════════════════
# PANEL 3 — ADMET detail
# ═══════════════════════════════════════════════════════════

admet <- data.frame(
  criterion = factor(
    c("PAINS exclusion","GI High absorption","LogS ≥ –6",
      "Lip. violations ≤ 1","All ADMET passed","BBB permeable\n(of 1,292)"),
    levels = rev(c("PAINS exclusion","GI High absorption","LogS ≥ –6",
                   "Lip. violations ≤ 1","All ADMET passed","BBB permeable\n(of 1,292)"))
  ),
  passed   = c(1363, 1464, 1387, 1464, 1292, 1136),
  total    = c(1464, 1464, 1464, 1464, 1464, 1292),
  is_highlight = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
)
admet$pct_pass <- round(admet$passed / admet$total * 100, 1)
admet$excl_label <- ifelse(
  admet$criterion == "PAINS exclusion",
  "101 PAINS excluded",
  paste0(admet$total - admet$passed, " excluded")
)

p3 <- ggplot(admet, aes(y = criterion)) +
  geom_col(aes(x = total), fill = "#EEEEEE", width = 0.65) +
  geom_col(aes(x = passed,
               fill = ifelse(is_highlight, col_cpd,
                      ifelse(criterion == "BBB permeable\n(of 1,292)",
                             "#5DCAA5", "#EF9F27"))),
           width = 0.65, show.legend = FALSE) +
  geom_text(aes(x = passed + 20,
                label = paste0(scales::comma(passed),
                               " (", pct_pass, "%)")),
            hjust = 0, size = 3.1, color = col_text, fontface = "bold") +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, 1750),
    breaks = c(0, 400, 800, 1200, 1600),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.22))
  ) +
  labs(
    x     = "Compounds passing",
    y     = NULL,
    title = "ADMET filter — per-criterion detail"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    plot.title         = element_text(size = 11, face = "bold", color = col_text,
                                      margin = margin(b = 6)),
    axis.text.y        = element_text(size = 9.5, color = col_text),
    axis.text.x        = element_text(size = 8.5, color = "#555555"),
    axis.title.x       = element_text(size = 8.5, color = "#555555"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.4),
    plot.background    = element_rect(fill = col_bg, color = NA),
    panel.background   = element_rect(fill = col_bg, color = NA)
  )

# ═══════════════════════════════════════════════════════════
# COMPOSE with patchwork
# ═══════════════════════════════════════════════════════════

layout <- "
AABB
AACC
"

combined <- p1 + p2 + p3 +
  plot_layout(design = layout) +
  plot_annotation(
    caption  = paste0(
      "Tukey threshold = Q1 + 1.5×IQR on final_score. ",
      "ADMET computed via RDKit (Python 3.10). ",
      "ChEMBL similarity ≥ 60%, pChEMBL ≥ 4.0."
    ),
    theme = theme(
      plot.caption    = element_text(size = 8, color = "#777777",
                                     hjust = 0, margin = margin(t = 8)),
      plot.background = element_rect(fill = col_bg, color = NA)
    )
  )

ggsave(
  "funnel_report_publication.png",
  plot   = combined,
  width  = 16,
  height = 14,
  dpi    = 300,
  bg     = col_bg
)

cat("Selesai: funnel_report_publication.png\n")
