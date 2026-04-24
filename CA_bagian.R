# ============================================================
# Correspondence Analysis — Plant Parts × Jamu Claim Categories
# Figure 4 in Tambunan et al. (2025)
# ============================================================
# Run from project root: Rscript scripts/02_visualization/CA_bagian.R

library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ── Baca data ────────────────────────────────────────────────
bagian <- read.csv("data/processed/jamu_ingredients_resolved.csv",
                   stringsAsFactors = FALSE)

# Standardize plant part names
bagian <- bagian %>%
  mutate(part_std = case_when(
    grepl("rhizom|rimpang",  bagian_tumbuhan, ignore.case=TRUE) ~ "Rhizome",
    grepl("daun|leaf|leave", bagian_tumbuhan, ignore.case=TRUE) ~ "Leaf",
    grepl("buah|fruit",      bagian_tumbuhan, ignore.case=TRUE) ~ "Fruit",
    grepl("kulit|bark",      bagian_tumbuhan, ignore.case=TRUE) ~ "Bark",
    grepl("biji|seed",       bagian_tumbuhan, ignore.case=TRUE) ~ "Seed",
    grepl("bunga|flower",    bagian_tumbuhan, ignore.case=TRUE) ~ "Flower",
    grepl("akar|root",       bagian_tumbuhan, ignore.case=TRUE) ~ "Root",
    grepl("seluruh|whole",   bagian_tumbuhan, ignore.case=TRUE) ~ "Whole plant",
    grepl("batang|stem",     bagian_tumbuhan, ignore.case=TRUE) ~ "Stem",
    TRUE ~ "Other"
  ))

# Build contingency matrix (plant part × claim category)
mat_parts <- bagian %>%
  filter(part_std != "Other") %>%
  group_by(part_std, category) %>%
  summarise(n = n_distinct(tanaman_accepted), .groups="drop") %>%
  pivot_wider(names_from=category, values_from=n, values_fill=0) %>%
  column_to_rownames("part_std") %>%
  as.matrix()

ca_parts <- CA(mat_parts, graph=FALSE)

# ── Plot ──────────────────────────────────────────────────────
row_df <- data.frame(
  name = rownames(ca_parts$row$coord),
  x    = ca_parts$row$coord[,1],
  y    = ca_parts$row$coord[,2]
)
col_df <- data.frame(
  name = rownames(ca_parts$col$coord),
  x    = ca_parts$col$coord[,1],
  y    = ca_parts$col$coord[,2]
)

eig  <- ca_parts$eig
pct1 <- round(eig[1,"percentage of variance"],1)
pct2 <- round(eig[2,"percentage of variance"],1)

p <- ggplot() +
  geom_hline(yintercept=0,color="grey70",linewidth=0.4,linetype="dashed") +
  geom_vline(xintercept=0,color="grey70",linewidth=0.4,linetype="dashed") +
  geom_point(data=row_df, aes(x=x,y=y),
             shape=21, fill="#2D6A4F", color="white",
             size=5, stroke=0.8) +
  geom_text_repel(data=row_df, aes(x=x,y=y,label=name),
                  size=3.5, fontface="bold", color="#1B4332",
                  point.padding=0.5, box.padding=0.6) +
  geom_segment(data=col_df, aes(x=0,y=0,xend=x*0.8,yend=y*0.8),
               arrow=arrow(length=unit(0.25,"cm"),type="closed"),
               color="#993556", linewidth=0.9) +
  geom_text(data=col_df, aes(x=x*0.9,y=y*0.9,label=name),
            size=4, fontface="bold", color="#993556") +
  labs(
    x     = sprintf("Dimension 1 (%.1f%%)", pct1),
    y     = sprintf("Dimension 2 (%.1f%%)", pct2),
    title = "Correspondence Analysis: Plant Parts × Jamu Claim Category",
    subtitle = sprintf("total inertia = %.1f%%", pct1+pct2)
  ) +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",color=NA))

ggsave("figures/fig4_CA_bagian.png", plot=p, width=10, height=9, dpi=300, bg="white")
cat("Saved: figures/fig4_CA_bagian.png\n")
