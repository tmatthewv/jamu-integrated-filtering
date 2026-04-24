# ============================================================
# Compound – Plant – Family Network  (versi bersih)
# Main graph: compound → tanaman → famili (yg punya senyawa)
# Non-priority plants: scatter kecil abu di bawah
# ============================================================

library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)

# ── 1. Baca data ─────────────────────────────────────────────
cpd_raw <- read.csv("compounds_final.csv",     stringsAsFactors = FALSE)
plants   <- read.csv("bobot_tanaman.csv",      stringsAsFactors = FALSE)
targets  <- read.csv("target_predictions.csv", stringsAsFactors = FALSE)

if ("SPECIES" %in% names(plants))
  names(plants)[names(plants) == "SPECIES"] <- "TANAMAN"

# ── 2. Compounds ─────────────────────────────────────────────
cpd <- cpd_raw %>%
  filter(compound_name != "Se-Methylselenocysteine") %>%
  select(compound_name, species, best_pchembl) %>%
  distinct(compound_name, species, .keep_all = TRUE)

tgt_map <- targets %>%
  filter(compound_name != "Se-Methylselenocysteine") %>%
  mutate(tgt_short = sub(".*\\((.+?)\\).*", "\\1", protein_name)) %>%
  group_by(compound_name) %>%
  summarise(target_label = paste(sort(unique(tgt_short)), collapse = "/"),
            .groups = "drop")

cpd <- cpd %>%
  left_join(tgt_map, by = "compound_name") %>%
  mutate(target_label = ifelse(is.na(target_label), "None", target_label))

tgt_pal <- c(AR="#4E79A7", ESR1="#D41E1E", HIF1A="#217A35",
             NFKB1="#8B3DAF", PTGS2="#E07B00", None="#AAAAAA")

cpd <- cpd %>%
  mutate(
    first_tgt = sub("/.*", "", target_label),
    cpd_color = tgt_pal[first_tgt],
    cpd_color = ifelse(is.na(cpd_color), "#AAAAAA", cpd_color),
    cpd_size  = 3.5 + pmin((best_pchembl - 4) * 1.5, 5)
  )

# ── 3. Priority plants ────────────────────────────────────────
Q1 <- quantile(plants$final_score, 0.25)
Q3 <- quantile(plants$final_score, 0.75)
tukey <- Q1 + 1.5 * (Q3 - Q1)
cat(sprintf("Tukey threshold : %.1f\n", tukey))

priority <- plants %>%
  filter(final_score >= tukey) %>%
  arrange(desc(final_score)) %>%
  group_by(FAMILI) %>% slice_head(n = 3) %>% ungroup()

non_priority <- plants %>%
  filter(!TANAMAN %in% priority$TANAMAN)

cat(sprintf("Priority plants : %d\n", nrow(priority)))
cat(sprintf("Non-priority    : %d\n", nrow(non_priority)))

# ── 4. Edges: HANYA nodes yang punya senyawa ─────────────────
e_cpd_plant <- cpd %>%
  filter(species %in% priority$TANAMAN) %>%
  transmute(from = compound_name, to = species, etype = "cpd_plant")

e_plant_fam <- priority %>%
  filter(TANAMAN %in% cpd$species) %>%   # hanya tanaman yg punya cpd
  transmute(from = TANAMAN, to = FAMILI, etype = "plant_fam")

edges_main <- bind_rows(e_cpd_plant, e_plant_fam) %>% distinct()

# ── 5. Nodes MAIN GRAPH ──────────────────────────────────────
cpd_nodes <- cpd %>%
  group_by(compound_name) %>% slice_head(n=1) %>% ungroup() %>%
  transmute(name=compound_name, ntype="Compound",
            nfill=cpd_color, nsize=cpd_size,
            label=compound_name, final_score=NA_real_)

plant_nodes <- priority %>%
  filter(TANAMAN %in% cpd$species) %>%
  transmute(name=TANAMAN, ntype="Plant_active",
            nfill="#2D6A4F",
            nsize=3.5 + pmin(final_score * 0.13, 5),
            label=sub("(\\S+) (\\S).*","\\1 \\2.", TANAMAN),
            final_score=final_score)

fam_nodes <- priority %>%
  filter(TANAMAN %in% cpd$species) %>%
  group_by(FAMILI) %>%
  summarise(fam_score=sum(final_score), .groups="drop") %>%
  transmute(name=FAMILI, ntype="Family",
            nfill="#F4A261",
            nsize=7 + pmin(fam_score * 0.04, 5),
            label=FAMILI, final_score=fam_score)

node_main <- bind_rows(cpd_nodes, plant_nodes, fam_nodes) %>%
  filter(name %in% unique(c(edges_main$from, edges_main$to)))

# ── 6. Priority plants tanpa senyawa (node hijau kecil) ──────
priority_no_cpd <- priority %>%
  filter(!TANAMAN %in% cpd$species)
cat(sprintf("Priority tanpa senyawa: %d\n", nrow(priority_no_cpd)))

# ── 7. Build graph ───────────────────────────────────────────
g <- graph_from_data_frame(edges_main, vertices=node_main, directed=FALSE)
E(g)$ecol   <- ifelse(edges_main$etype=="cpd_plant", "#444444", "#999999")
E(g)$ealpha <- ifelse(edges_main$etype=="cpd_plant", 0.75, 0.45)
E(g)$ewidth <- ifelse(edges_main$etype=="cpd_plant", 1.0,  0.6)

set.seed(42)
lay <- create_layout(g, layout="stress")

xr <- range(lay$x); yr <- range(lay$y)
cat(sprintf("Layout range: x=[%.2f,%.2f] y=[%.2f,%.2f]\n",
            xr[1],xr[2],yr[1],yr[2]))

# ── 8. Non-priority coords (grid di bawah) ───────────────────
n_np    <- nrow(non_priority)
nc_np   <- 10
nr_np   <- ceiling(n_np / nc_np)
xstep   <- (xr[2]-xr[1]) / (nc_np - 1)
ybase   <- yr[1] - 2.0

np_df <- data.frame(
  name   = non_priority$TANAMAN,
  FAMILI = non_priority$FAMILI,
  col_i  = (seq_len(n_np)-1) %% nc_np,
  row_i  = floor((seq_len(n_np)-1) / nc_np)
) %>% mutate(
  x_plot = xr[1] + col_i * xstep,
  y_plot = ybase - row_i * 0.42
)

# Priority plants tanpa senyawa → row terpisah di bawah graph
n_pnc  <- nrow(priority_no_cpd)
if (n_pnc > 0) {
  pnc_x <- seq(xr[1], xr[2], length.out=n_pnc)
  priority_no_cpd <- priority_no_cpd %>%
    mutate(x_plot = pnc_x, y_plot = yr[1] - 1.0)
}

# ── 9. Legend data ────────────────────────────────────────────
leg_x <- xr[2] + 0.8
leg_y <- seq(yr[2], yr[2]-2.5, length.out=9)
leg_df <- data.frame(
  x     = leg_x,
  y     = leg_y,
  label = c("Compound target:","  AR","  ESR1","  HIF1A",
            "  NF-kB1 (NFKB1)","  COX-2 (PTGS2)",
            "","Priority plant","Family"),
  color = c(NA,"#4E79A7","#D41E1E","#217A35","#8B3DAF","#E07B00",
            NA,"#2D6A4F","#F4A261"),
  shape = c(NA, 21,21,21,21,21, NA, 21, 23),
  sz    = c(NA, 3,3,3,3,3, NA, 4,4),
  is_pt = c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE,TRUE,TRUE)
)

# ── 10. Plot ─────────────────────────────────────────────────
p <- ggraph(lay) +

  geom_edge_link(
    aes(colour=ecol, alpha=ealpha, width=ewidth),
    show.legend=FALSE
  ) +
  scale_edge_colour_identity() +
  scale_edge_alpha_identity() +
  scale_edge_width_identity() +

  # Non-priority plants (scatter abu)
  geom_point(data=np_df, aes(x=x_plot, y=y_plot),
             shape=21, fill="#CCCCCC", color="#AAAAAA",
             size=1.8, stroke=0.25, alpha=0.5,
             inherit.aes=FALSE) +

  # Priority tanpa senyawa (hijau, ukuran sedang)
  { if (n_pnc > 0)
      geom_point(data=priority_no_cpd,
                 aes(x=x_plot, y=y_plot),
                 shape=21, fill="#52B788", color="#2D6A4F",
                 size=3.5, stroke=0.6, alpha=0.7,
                 inherit.aes=FALSE)
  } +

  # Family (diamond)
  geom_node_point(
    data=function(x) x[x$ntype=="Family",],
    aes(size=nsize),
    shape=23, fill="#F4A261", color="#C65C1A", stroke=0.9
  ) +

  # Priority plants dengan senyawa
  geom_node_point(
    data=function(x) x[x$ntype=="Plant_active",],
    aes(size=nsize),
    shape=21, fill="#2D6A4F", color="#1B4332", stroke=0.8
  ) +

  # Compounds
  geom_node_point(
    data=function(x) x[x$ntype=="Compound",],
    aes(size=nsize, fill=nfill),
    shape=21, color="white", stroke=0.9
  ) +
  scale_fill_identity() +
  scale_size_identity() +

  # Label Family
  geom_node_text(
    data=function(x) x[x$ntype=="Family",], aes(label=label),
    fontface="bold", size=3.5, color="#7B3000",
    repel=TRUE, max.overlaps=30,
    point.padding=0.6, box.padding=0.7,
    segment.color="#C65C1A", segment.size=0.35,
    show.legend=FALSE
  ) +

  # Label Plants
  geom_node_text(
    data=function(x) x[x$ntype=="Plant_active",], aes(label=label),
    fontface="italic", size=2.9, color="#1B4332",
    repel=TRUE, max.overlaps=30,
    point.padding=0.5, box.padding=0.6,
    segment.color="#2D6A4F", segment.size=0.3,
    show.legend=FALSE
  ) +

  # Label Compounds
  geom_node_text(
    data=function(x) x[x$ntype=="Compound",], aes(label=label),
    size=2.7, color="#111111",
    repel=TRUE, max.overlaps=40,
    point.padding=0.5, box.padding=0.6,
    segment.color="#666666", segment.size=0.25,
    show.legend=FALSE
  ) +

  # ── Legend ───────────────────────────────────────────────
  annotate("text", x=leg_x, y=yr[2]+0.15,
           label="Legend", hjust=0, size=3.4,
           fontface="bold", color="#222222") +
  annotate("point",
           x=leg_df$x[leg_df$is_pt],
           y=leg_df$y[leg_df$is_pt],
           shape=leg_df$shape[leg_df$is_pt],
           fill=leg_df$color[leg_df$is_pt],
           color=ifelse(leg_df$shape[leg_df$is_pt]==23,"#C65C1A","white"),
           size=leg_df$sz[leg_df$is_pt], stroke=0.7) +
  annotate("text",
           x=leg_df$x + 0.18, y=leg_df$y,
           label=leg_df$label,
           hjust=0, size=2.9, color="#333333") +

  # Non-priority label
  annotate("text", x=xr[1], y=ybase+0.3,
           label=sprintf("Non-priority plants  (n=%d, final_score < %.0f)",
                         n_np, tukey),
           hjust=0, size=2.7, color="#888888", fontface="italic") +
  annotate("point", x=xr[1]-0.2, y=ybase+0.3,
           shape=21, fill="#CCCCCC", color="#AAAAAA",
           size=2.0, stroke=0.3) +

  theme_void(base_size=11) +
  theme(
    plot.title      = element_text(size=14, face="bold", hjust=0.5,
                                   margin=margin(b=4)),
    plot.subtitle   = element_text(size=9, hjust=0.5, color="grey40",
                                   margin=margin(b=8)),
    plot.background = element_rect(fill="white", color=NA),
    plot.margin     = margin(20, 110, 20, 20)
  ) +
  labs(
    title    = "Compound \u2013 Plant \u2013 Family Network",
    subtitle = sprintf(
      "%d bioactive compounds  \u00b7  %d priority plants  \u00b7  %d families  \u00b7  %d non-priority plants",
      n_distinct(cpd$compound_name), nrow(priority),
      n_distinct(priority$FAMILI), nrow(non_priority)
    )
  )

ggsave("network_compound_plant_family.png",
       plot=p, width=20, height=16, dpi=300, bg="white")
cat("Saved: network_compound_plant_family.png\n")
