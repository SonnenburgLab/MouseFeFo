# Metabolomics Analysis
#
# Figure assignment:
#   1A - Schematic
#   1B - C2 vs fermentation time (triplicates averaged BEFORE PPCA)
#   1C - Volcano: final vs baseline, pooled across all vegetable foods
#   1D - Heatmap of log2FC (final / baseline) per food x metabolite

#Prepare tables####
rm(list = ls())

source("~/Documents/FeFoPaper_Analysis/0-config.R")

# Read in and load data #####
mapping_file <- read.csv("02082023_tcFeFo_mapping_dupl.csv",
  stringsAsFactors = FALSE) %>%
  mutate(Timepoint_d = as.numeric(Timepoint_d),
         Timepoint_h = as.numeric(Timepoint_h))

clean_metab <- read.csv("20250221_clean_metabolites.csv",
  stringsAsFactors = FALSE) 

metab_class <- read.csv("04202026_cleanmetab_class.csv")

metab_class <- metab_class %>%
  mutate(Superclass = ifelse(is.na(Superclass), "Other / Unknown", Superclass))

# Vegetable sample selection #####
pca_samples <- mapping_file %>%
  filter(Category == "Vegetable",
         Type     != "Starter") %>%
  pull(Sample_name) %>% base::unique()

pca_long <- clean_metab %>%
  filter(Sample_name %in% pca_samples) %>% 
  filter(library == "Sonn") %>% 
  distinct(Sample_name, Metabolite, Metabolite_mode, .keep_all = TRUE)

# FIGURE 1B — PC2 vs Fermentation Time #####
## STEP 1 — Average triplicates BEFORE PPCA####
avg_long <- pca_long %>%
  mutate(group = sub("_trip[0-9]+$", "", Sample_name)) %>%
  group_by(group, Metabolite_mode) %>%
  summarise(peak_area_avg = mean(peak_area_norm, na.rm = TRUE), .groups = "drop")

avg_wide <- avg_long %>%
  pivot_wider(names_from  = Metabolite_mode,
              values_from = peak_area_avg) #%>% 
  #select(-c(" X_c18pos"))

avg_ids  <- avg_wide$group
pca_mat  <- as.matrix(avg_wide[, -1])
rownames(pca_mat) <- avg_ids

min_nonzero <- min(pca_mat[pca_mat > 0], na.rm = TRUE)
pca_mat[is.na(pca_mat) | pca_mat == 0] <- min_nonzero / 2
pca_mat <- log2(pca_mat)
pca_mat <- pca_mat[, apply(pca_mat, 2, var) > 0]

message("PPCA matrix (averaged): ", nrow(pca_mat), " x ", ncol(pca_mat))

## STEP 2 — Probabilistic PCA ####
nPcs_use <- min(10, nrow(pca_mat) - 1)
pca_res  <- pcaMethods::pca(pca_mat, method = "ppca",
                            nPcs = nPcs_use, scale = "uv", center = TRUE)
var_exp  <- round(pca_res@R2 * 100, 1)

message("Variance explained: ", paste0(var_exp, "%", collapse = ", "))

ylab_pc1 <- paste0("PC1 (", var_exp[1], "% variance explained)")
ylab_pc2 <- paste0("PC2 (", var_exp[2], "% variance explained)")

pca_scores <- as.data.frame(pcaMethods::scores(pca_res))
pca_scores$group <- rownames(pca_scores)

time_data <- mapping_file %>%
  filter(Category == "Vegetable",
         Type     != "Starter",
         !Sample_name %in% "Curtido_F_T96_trip3") %>%
  mutate(group = sub("_trip[0-9]+$", "", Sample_name)) %>%
  distinct(group, Label, .keep_all = TRUE) %>%
  inner_join(pca_scores, by = "group") %>%
  arrange(Label, Timepoint_d) %>%
  mutate(Timepoint_d_fac = factor(as.character(Timepoint_d),
                                  levels = as.character(sort(base::unique(Timepoint_d)))))

feat_label <- paste0(n_distinct(pca_long$Metabolite_mode), " annotated features")

## STEP 3 — Make Plot ####
p_1B <- ggplot(time_data,
               aes(x = Timepoint_d_fac, y = PC2, color = Label)) +
  geom_line(aes(group = Label), alpha = 0.8, linewidth = 2) +
  geom_point(size = 5.5, alpha = 0.9) +
  scale_color_manual(values = vegetable_colors, na.value = "grey60") +
  labs(
    title    = "",
    x        = "Fermentation (d)",
    y        = ylab_pc2,
    color    = NULL
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 11) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    axis.text        = element_text(color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4)
  )

# FIGURE 1C — Volcano Plot (final vs baseline, all vegetable foods) #####
results <- clean_metab %>% 
  filter(library == "Sonn") %>% 
  left_join(mapping_file, relationship = "many-to-many") %>% 
  filter(!is.na(Order)) %>% 
  group_by(Label, Metabolite_mode, library) %>%
  summarise(
    mean_F = mean(peak_area_norm[Order == "F"], na.rm = TRUE),
    mean_I = mean(peak_area_norm[Order == "I"], na.rm = TRUE),
    difference = mean_F - mean_I,
    fold_change = mean_F / mean_I,
    log2FC = log2(mean_F / mean_I),
    p_value = tryCatch(
      t.test(peak_area_norm[Order == "F"],
             peak_area_norm[Order == "I"])$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%  group_by(Label) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup()

top_label <- volcano_stats %>%
  filter(Metabolite_mode %in% c(
    "2-HYDROXY-3-METHYLPENTANOIC ACID_c18neg",
    "3-PHENYLLACTIC ACIDPEAK2_c18neg",
    "HYDROXYPHENYLLACTIC ACID_c18neg",
    "LEUCINE_c18pos",
    "PHENYLALANINEPEAK1_c18neg",
    "STACHYOSE_c18neg",
    "TRYPTOPHAN_c18pos",
    "TRYPTAMINE_c18pos",
    "MALIC ACID_c18pos",
    "ADENOSINE_c18pos",
    "INOSINE_c18pos"
  )) %>%
  group_by(Metabolite) %>%
  slice_max(neg_log10p, n = 4, with_ties = FALSE) %>%
  ungroup()

# Merge with metabolite class info
volcano_stats <- results %>%
  mutate(Metabolite = sub("_c18neg$|_c18pos$", "", Metabolite_mode)) %>%
  left_join(metab_class %>% distinct(Metabolite, Superclass), by = "Metabolite") %>%
  mutate(
    neg_log10p = -log10(p_adj),
    significant = p_adj < 0.05 & abs(log2FC) >= 2
  )

p_volcano <- ggplot(volcano_stats,
                    aes(x = log2FC, y = neg_log10p)) +
  geom_point(data = \(d) filter(d, !significant),
             color = "grey60", alpha = 0.4, size = 1.8) +
  geom_point(data = \(d) filter(d, significant),
             aes(color = Superclass), alpha = 1.0, size = 3) +
  scale_color_manual(values = metabolite_class_colors, na.value = "grey60",
                     name = "Metabolite class") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = top_label,
    aes(label = Metabolite),
    size = 1.5, max.overlaps = 30, segment.color = "grey60",
    show.legend = FALSE
  ) +
  labs(
    title    = "Metabolite changes: final vs baseline",
    subtitle = paste0("Sonn library  |  BH-adj p < 0.05  |  |log2FC| \u2265 2  |  top 20 labeled"),
    x        = "log2 fold change (final / baseline)",
    y        = "-log10(p-adj)"
  ) +
  theme_minimal(base_family = "Helvetica", base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position  = "right",
    axis.text        = element_text(color = "black"),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4)
  )

# FIGURE 1D — Heatmap: log2FC per food x metabolite (all Sonn library metabolites) #####
library(ComplexHeatmap)
library(circlize)

# Build matrix from results, Sonn only
lfc_t <- results %>%
  filter(library == "Sonn") %>%
  mutate(Metabolite = sub("_c18neg$|_c18pos$", "", Metabolite_mode)) %>%
  select(Label, Metabolite, log2FC) %>%
  pivot_wider(names_from = Metabolite, values_from = log2FC) %>%
  column_to_rownames("Label") %>%
  as.matrix()

lfc_t_scaled <- scale(lfc_t)

# Color scale
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2),
                      c("#0072b2", "#56b4e9", "lightyellow", "#E69F00", "#d55e00"))

# Left annotation: bar chart of n metabolites with z-score > 1.5 per food
n_elevated_per_food <- rowSums(lfc_t_scaled > 1, na.rm = TRUE)

row_ann_1D <- ComplexHeatmap::rowAnnotation(
  `z > 1.5` = ComplexHeatmap::anno_barplot(
    n_elevated_per_food,
    bar_width = 0.8,
    gp        = grid::gpar(
      fill = ifelse(names(n_elevated_per_food) == "Green Sauerkraut (wNS)",
                    "#99c13d", "grey75")
    ),
    width      = grid::unit(2, "cm"),
    axis_param = list(gp = grid::gpar(fontsize = 7))
  )
)

# Heatmap
heatmap_bar <- ComplexHeatmap::Heatmap(
  lfc_t_scaled,
  col                         = col_fun,
  name                        = "log2FC",
  left_annotation             = row_ann_1D,
  cluster_rows                = TRUE,
  cluster_columns             = TRUE,
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  row_km                      = 3,
  show_row_names              = TRUE,
  show_column_names           = TRUE,
  row_names_gp                = grid::gpar(fontsize = 9),
  column_names_gp             = grid::gpar(fontsize = 7),
  na_col                      = "grey90",
  rect_gp                     = grid::gpar(col = NA),
  column_title                = paste0("log2FC (final / baseline)  |  ",
                                       ncol(lfc_t), " metabolites")
)