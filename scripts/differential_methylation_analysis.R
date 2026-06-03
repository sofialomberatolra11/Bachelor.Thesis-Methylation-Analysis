
### DIFFERENTIAL METHYLATION ANALYSIS ###

# Required packages
library(limma)
library(dplyr)
library(pheatmap)
library(grid)
library(IlluminaHumanMethylationEPICv2anno.ilm10b4.hg38)

### LOAD DATA ###

load("data/DNAm_matrix.RData")
load("data/group_df.RData")

### MATCH SAMPLE ORDER ###

matriu.beta <- matriu.beta[
  ,
  match(
    group_df$Sample,
    colnames(matriu.beta)
  )
]

### DEFINE GROUPS ###

group_df$Group <- as.factor(group_df$Group)

# Create design matrix
model.s <- model.matrix(
  ~0 + Group,
  data = group_df
)

colnames(model.s) <- paste0(
  "cluster",
  levels(group_df$Group)
)

### FIT LINEAR MODEL ###

fit.s <- lmFit(
  matriu.beta,
  model.s
)

### CREATE CONTRAST MATRIX ###

contMatrix.s <- makeContrasts(
  cluster1 - cluster2,
  levels = model.s
)

### APPLY CONTRASTS ###

fit2.s <- contrasts.fit(
  fit.s,
  contMatrix.s
)

fit2.s <- eBayes(fit2.s)

### LOAD EPIC ARRAY ANNOTATION ###

annEPICv2 <- getAnnotation(
  IlluminaHumanMethylationEPICv2anno.20a1.hg38
)

### EXTRACT DIFFERENTIAL CpGs ###

diffClusters.s <- topTable(
  fit = fit2.s,
  num = Inf,
  coef = 1
)

# Add methylation difference
diffClusters.s$deltaB <- fit2.s$coefficients[
  match(
    rownames(diffClusters.s),
    rownames(fit2.s$coefficients)
  )
]

### FILTER SIGNIFICANT CpGs ###

SigDiffClusters <- diffClusters.s %>%
  filter(
    adj.P.Val < 0.01,
    abs(deltaB) > 0.25
  )

# Extract probe names
differential_probes <- rownames(
  SigDiffClusters
)

differential_probes <- sub(
  "_.*",
  "",
  differential_probes
)

### FILTER NORMAL SAMPLES ###

normals_epic.filter <- normals_epic[
  rownames(normals_epic) %in%
    differential_probes,
]

### FILTER PBL SAMPLES ###

rownames(matriu.beta) <- sub(
  "_.*",
  "",
  rownames(matriu.beta)
)

matriu.beta.filtered <- matriu.beta[
  rownames(matriu.beta) %in%
    differential_probes,
]

### FIND COMMON PROBES ###

common_probes <- intersect(
  rownames(normals_epic.filter),
  rownames(matriu.beta.filtered)
)

# Keep common probes
normals_epic.filter <- normals_epic.filter[
  common_probes,
]

matriu.beta.filtered <- matriu.beta.filtered[
  common_probes,
]

### CONVERT TO MATRICES ###

tumor_matrix <- as.matrix(
  matriu.beta.filtered
)

normal_matrix <- as.matrix(
  normals_epic.filter
)

### HEATMAP - PBL SAMPLES ###

hm_tumor <- pheatmap(
  tumor_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  use_raster = FALSE
)

### DISPLAY HEATMAP ###

grid.newpage()
grid.draw(hm_tumor$gtable)

### HEATMAP - NORMAL B CELLS ###

hm_normal <- pheatmap(
  normal_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  use_raster = FALSE
)

### ORDER NORMAL CELL TYPES ###

correct_order <- c(
  "NBC.NC11_42_NBC",
  "NBC.N11_54_NBC",
  "GC.T15_26_GC",
  "MBC.NC11_42_MBC",
  "MBC.N11_80_ncsMBC",
  "MBC.NC11_55_ncsMBC",
  "PC.T15_26_PC"
)

normal_matrix <- normal_matrix[
  ,
  correct_order
]

### HEATMAP - ORDERED NORMAL CELLS ###

hm_normal_ordered <- pheatmap(
  normal_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100)
)

### COMBINED HEATMAP ###

combined_matrix <- cbind(
  normal_matrix,
  tumor_matrix
)

hm_combined <- pheatmap(
  combined_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(
    c("blue", "white", "red")
  )(100),
  use_raster = FALSE
)

