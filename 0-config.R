# =============================================================================
# 0-config.R — Project-Wide Configuration
# FeFoPaper Analysis
# -----------------------------------------------------------------------------
# SOURCE THIS FILE at the top of every figure script:
#   source("0-config.R")
# =============================================================================


# --- Working Directory --------------------------------------------------------

# --- Packages -----------------------------------------------------------------
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(rstatix)
library(patchwork)
library(svglite)
library(DESeq2)
library(ggrepel)
library(conflicted)
library(pcaMethods)
library(pheatmap)
library(ggrepel)
library(pcaMethods)
library(pheatmap)
library(colorRamp2)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2",      quietly = TRUE)) BiocManager::install("DESeq2")
if (!requireNamespace("ggrepel",     quietly = TRUE)) install.packages("ggrepel")
conflict_prefer("select",    "dplyr")
conflict_prefer("filter",    "dplyr")
conflict_prefer("rename",    "dplyr")
conflict_prefer("intersect", "base")
conflict_prefer("union",     "base")
conflict_prefer("setdiff",   "base")
conflict_prefer("unname",    "base")
conflict_prefer("%in%", "base")
conflict_prefer("unique", "base")
conflicts_prefer(dplyr::slice)


# Auto-install any packages not yet installed
installed <- rownames(installed.packages())
to_install <- required_packages[!required_packages %in% installed]
if (length(to_install) > 0) {
  message("Installing missing packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))


# --- Color Palette ------------------------------------------------------------
# Primary diet group colors — use these in every figure
diet_colors <- c(
  "Water"    = "#cee0ee",   # light blue  — water/control group
  "SK-FDMs"  = "#84a45f"    # olive green — 3-week fermented diet group
)
# Alias so either name works in scripts
group_colors <- diet_colors

# Figure 3 microbiome scripts use the label "3wSK" rather than "SK-FDMs".
# These aliases expose the same hex values under the names those scripts expect,
# so no script needs to redefine colours locally.
col_water     <- "#cee0ee"   # water / control
col_3wSK      <- "#84a45f"   # 3-week fermented diet
col_unferment <- "#d4956a"   # HFD Unfermented group (Figure 3 HFD only)

diet_cols     <- c("Water" = col_water, "3wSK"         = col_3wSK)

# GF takedown physiology figures — slightly lighter palette for germ-free mice
gf_col_water  <- "#e6eff6"   # light blue  (GF water/control)
gf_col_3wSK   <- "#c1d1af"   # olive green (GF 3-week fermented)
gf_diet_cols  <- c("Water" = gf_col_water, "3wSK" = gf_col_3wSK)
diet_cols_hfd <- c("Water" = col_water,
                   "Unfermented"  = col_unferment,
                   "3w Fermented" = col_3wSK)

# Vegetable fermentation colors — keyed by Label column in tcFeFo mapping
# Used in Figure 1 metabolomics plots
vegetable_colors <- c(
  "Baechu Kimchi (wNS)"    = "#cd5c5c",
  "Baechu Kimchi (wNVS)"   = "#f2673a",
  "Baechu Kimchi (wS)"     = "#6d0202",
  "Baechu Kimchi (wVS)"    = "#d73017",
  "Curtido"                = "#d994a6",
  "Green Sauerkraut (wNS)" = "#99c13d",
  "Green Sauerkraut (wS)"  = "#679e33",
  "Nabak Kimchi (wNS)"     = "#c09510",
  "Nabak Kimchi (wS)"      = "#dbc06d",
  "Red Sauerkraut"         = "#71265d"
)

# Fermented food type colors — keyed by Ferment column in tcFeFo mapping
# Used in Figure 1 metabolomics PCA plots
ferment_colors <- c(
  "Kimchi"    = "#d73017",   # red
  "Curtido"   = "#d994a6",   # pink
  "Sauekraut" = "#679e33"    # green (matches spelling in mapping file)
)

# Metabolite Superclass colors — ClassyFire hierarchy
# Used in Figure 1C volcano plot
metabolite_class_colors <- c(
  "Organic acids and derivatives"              = "#e07b54",
  "Lipids and lipid-like molecules"            = "#d19779",
  "Organoheterocyclic compounds"               = "#719bae",
  "Nucleosides, nucleotides, and analogues"    = "#f5b27f",
  "Benzenoids"                                 = "#ba883c",
  "Phenylpropanoids and polyketides"           = "#4f7d9b",
  "Organic oxygen compounds"                   = "#2a8787",
  "Organic nitrogen compounds"                 = "#34495e",
  "Alkaloids and derivatives"                  = "#e67e22",
  "Carbohydrates and carbohydrate conjugates"  = "#f39c12",
  "Other / Unknown"                            = "#3d5e7d"
)

col_fun = colorRamp2(c(-2, -1, 0, 1, 2), c("#0072b2","#56b4e9", "lightyellow", "#E69F00", "#d55e00"))
col_fun(seq(-3, 3))

ferment_colors <- c(
  "Amasi" = "#99b9ff", 
  "Kimchi" = "#5f918d", 
  "Curtido" = "#8f925a", 
  "Filmjolk" = "#d2d7ee", 
  "Kefir" = "#3372A9", 
  "Sauerkraut" = "#8ab771", 
  "Viili" = "#8ac0e5", 
  "Yogurt" = "#73accb"
)

# --- ggplot2 Theme ------------------------------------------------------------
# Applied to all figures; override per-plot with + theme() as needed
base_theme <- theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y     = element_text(color = "black"),
    axis.line       = element_line(color = "black", linewidth = 0.5),
    axis.ticks      = element_line(color = "black", linewidth = 0.5),
    plot.title      = element_text(size = 11, face = "bold")
  )

# Apply as the global default so you don't need to add it every time
# (can still override per-plot with + theme(...))
theme_set(theme_classic(base_size = 11))


# --- Significance Notation ----------------------------------------------------
# Used in add_significance() and stat_pvalue_manual() calls
sig_cutpoints <- c(0, 0.001, 0.01, 0.05, 1)
sig_symbols   <- c("***", "**", "*", "n.s.")


# --- Figure Output Dimensions (Nature Communications) -------------------------
# Single-column panel:  89mm wide  ≈ 3.50 in
# Double-column figure: 183mm wide ≈ 7.20 in
# Typical panel height: 3.0–4.0 in depending on content
fig_width_single  <- 3.50   # inches
fig_width_double  <- 7.20   # inches
fig_height_default <- 3.50  # inches
fig_dpi            <- 300   # minimum for raster output


# --- Diet Label Recoding ------------------------------------------------------
# Canonical mapping from raw data values → display labels used in figures
# Add new entries here as new diet groups appear in the data
diet_labels <- c(
  "3w"    = "SK-FDMs",
  "Water" = "Water"
)

# Factor level order for consistent x-axis ordering across all plots
diet_levels <- c("Water", "SK-FDMs")
