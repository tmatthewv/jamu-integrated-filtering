# ============================================================
# Correspondence Analysis — Plant Species × Jamu Claim Categories
# Figure 5 in Tambunan et al. (2025)
# ============================================================
# Run from project root: Rscript scripts/02_visualization/CA_tanaman.R

library(FactoMineR)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ── Baca data ────────────────────────────────────────────────
plants <- read.csv("data/raw/bobot_tanaman.csv", stringsAsFactors = FALSE)
if ("SPECIES" %in% names(plants)) names(plants)[names(plants)=="SPECIES"] <- "TANAMAN"

# Filter: total_freq >= 5
active <- plants %>% filter((MALE_FREQ + FEM_FREQ + NEU_FREQ) >= 5)
cat(sprintf("Plants in CA: %d\n", nrow(active)))

# ── Contingency matrix ────────────────────────────────────────
mat <- as.matrix(active[, c("MALE_FREQ","FEM_FREQ","NEU_FREQ")])
rownames(mat) <- active$TANAMAN

# ── CA ────────────────────────────────────────────────────────
ca_res <- CA(mat, graph = FALSE)

# Export koordinat untuk network visualization
coords <- data.frame(
  TANAMAN = rownames(ca_res$row$coord),
  Dim1    = ca_res$row$coord[,1],
  Dim2    = ca_res$row$coord[,2]
)
write.csv(coords, "data/processed/tanaman_CA_koordinat.csv", row.names = FALSE)

# ── Plot ──────────────────────────────────────────────────────
row_df <- data.frame(
  name   = rownames(ca_res$row$coord),
  x      = ca_res$row$coord[,1],
  y      = ca_res$row$coord[,2],
  size   = active$MALE_FREQ + active$FEM_FREQ + active$NEU_FREQ,
  family = active$FAMILI
) %>% left_join(active %>% select(TANAMAN, FAMILI), by = c("name"="TANAMAN"))

col_df <- data.frame(
  name = rownames(ca_res$col$coord),
  x    = ca_res$col$coord[,1],
  y    = ca_res$col$coord[,2]
)

eig <- ca_res$eig
pct1 <- round(eig[1,"percentage of variance"], 1)
pct2 <- round(eig[2,"percentage of variance"], 1)

top_fam <- row_df %>% count(FAMILI.y, wt=size) %>%
  top_n(8, n) %>% pull(FAMILI.y)
row_df <- row_df %>%
  mutate(fam_label = ifelse(FAMILI.y %in% top_fam, FAMILI.y, "Other"))

fam_pal <- c(
  "Zingiberaceae"="#2D6A4F","Piperaceae"="#4E79A7","Apiaceae"="#F28E2B",
  "Myrtaceae"="#E15759","Amaryllidaceae"="#76B7B2","Lauraceae"="#59A14F",
  "Myristicaceae"="#EDC948","Asteraceae"="#B07AA1","Other"="#CCCCCC"
)

p <- ggplot() +
  geom_hline(yintercept=0, color="grey70", linewidth=0.4, linetype="dashed") +
  geom_vline(xintercept=0, color="grey70", linewidth=0.4, linetype="dashed") +
  geom_point(data=row_df, aes(x=x, y=y, fill=fam_label, size=size),
             shape=21, color="white", stroke=0.6, alpha=0.85) +
  geom_segment(data=col_df, aes(x=0,y=0,xend=x*0.8,yend=y*0.8),
               arrow=arrow(length=unit(0.25,"cm"),type="closed"),
               color="#993556", linewidth=0.9) +
  geom_text(data=col_df, aes(x=x*0.88, y=y*0.88, label=name),
            size=3.8, fontface="bold", color="#993556") +
  geom_text_repel(data=row_df %>% filter(size >= 10),
                  aes(x=x, y=y, label=sub("(\\S+) (\\S).*","\\1 \\2.",name)),
                  size=2.5, fontface="italic", max.overlaps=30,
                  point.padding=0.3, box.padding=0.5,
                  segment.color="grey50", segment.size=0.3) +
  scale_fill_manual(values=fam_pal, name="Family") +
  scale_size_continuous(range=c(2,9), name="Total freq") +
  labs(
    x       = sprintf("Dimension 1 (%.1f%%)", pct1),
    y       = sprintf("Dimension 2 (%.1f%%)", pct2),
    title   = "Correspondence Analysis: Plant Species × Jamu Claim Category",
    subtitle= sprintf("n=%d species  ·  set.seed(42)  ·  total inertia explained = %.1f%%",
                      nrow(row_df), pct1+pct2)
  ) +
  theme_minimal(base_size=11) +
  theme(plot.background=element_rect(fill="white",color=NA))

ggsave("figures/fig5_CA_tanaman.png", plot=p, width=12, height=10, dpi=300, bg="white")
cat("Saved: figures/fig5_CA_tanaman.png\n")
