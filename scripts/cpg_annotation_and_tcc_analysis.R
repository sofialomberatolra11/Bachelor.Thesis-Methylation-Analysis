
### CpG ANNOTATION AND TCC ANALYSIS ###

# Required packages
library(minfi)
library(minfiData)
library(dbplyr)
library(dplyr)
library(ggplot2)
library(viridis)

### FILTER RGSET SAMPLES ###

# Keep valid patient samples
RGSET <- RGSet.all[, colnames(RGSet.all) %in% cols]

# Remove duplicated samples
duplicated_samples <- c("PBL_76T", "PBL_72T", "PBL_108T")

RGSet.filtered <- RGSET[
  ,
  !colnames(RGSET) %in% duplicated_samples
]

### PREPARE PROBE ANNOTATION ###

# Assign correct CpG probe names
probe_names <- data$V1

# Preprocess raw methylation data
Mset.pre <- preprocessRaw(RGSet.filtered)

# Keep probes present in the dataset
Mset.filtered <- Mset.pre[
  rownames(Mset.pre) %in% probe_names,
]

# Retrieve annotation
annot <- getAnnotation(Mset.filtered)

### CREATE BETA MATRIX WITH CpG IDs ###

data.annotation <- data

cols <- grep(
  "^PBL|V1",
  colnames(data.annotation),
  value = TRUE
)

data.annotation <- data[, ..cols]

annot.beta <- annot[
  rownames(annot) %in% probe_names,
]

### MAP INTERNAL IDS TO REAL CpG PROBE IDS ###

mapping <- data.frame(
  internal_id = rownames(data),
  real_probe = data[, "V1"]
)

### ANNOTATION - SAMPLE SET 76 ###

filtered_probes_76 <- mapping$real_probe[
  mapping$internal_id %in% sondes_filtrades.76
]

annot_filtered_76 <- annot[
  rownames(annot) %in% filtered_probes_76,
]

info_probes_76 <- annot_filtered_76[
  ,
  c(
    "chr",
    "pos",
    "Relation_to_Island",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Group",
    "Islands_Name"
  )
]

head(info_probes_76)

### ANNOTATION - SAMPLE SET 72 ###

filtered_probes_72 <- mapping$real_probe[
  mapping$internal_id %in% sondes_filtrades.72
]

annot_filtered_72 <- annot[
  rownames(annot) %in% filtered_probes_72,
]

info_probes_72 <- annot_filtered_72[
  ,
  c(
    "chr",
    "pos",
    "Relation_to_Island",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Group",
    "Islands_Name"
  )
]

head(info_probes_72)

### ANNOTATION - SAMPLE SET 108 ###

filtered_probes_108 <- mapping$real_probe[
  mapping$internal_id %in% sondes_filtrades.108
]

annot_filtered_108 <- annot[
  rownames(annot) %in% filtered_probes_108,
]

info_probes_108 <- annot_filtered_108[
  ,
  c(
    "chr",
    "pos",
    "Relation_to_Island",
    "UCSC_RefGene_Name",
    "UCSC_RefGene_Group",
    "Islands_Name"
  )
]

head(info_probes_108)

### COMBINE PROBE SETS ###

vec_76 <- as.character(info_probes_76$probe)
vec_72 <- as.character(info_probes_72$probe)
vec_108 <- as.character(info_probes_108$probe)

combined_vec <- unique(
  c(vec_76, vec_72, vec_108)
)

### RELATION TO GENE REGIONS ###

annot_df <- as.data.frame(annot)

info72_df <- as.data.frame(info_probes_72)
info76_df <- as.data.frame(info_probes_76)
info108_df <- as.data.frame(info_probes_108)

df_all <- bind_rows(
  annot_df %>%
    select(UCSC_RefGene_Group) %>%
    mutate(group = "Background"),
  
  info72_df %>%
    select(UCSC_RefGene_Group) %>%
    mutate(group = "Set72"),
  
  info76_df %>%
    select(UCSC_RefGene_Group) %>%
    mutate(group = "Set76"),
  
  info108_df %>%
    select(UCSC_RefGene_Group) %>%
    mutate(group = "Set108")
)

# Define plotting order
df_all$group <- factor(
  df_all$group,
  levels = c(
    "Background",
    "Set72",
    "Set76",
    "Set108"
  )
)

### COUNT GENE REGION FREQUENCIES ###

df_plot_gene <- df_all %>%
  count(
    group,
    UCSC_RefGene_Group,
    name = "n"
  ) %>%
  rename(
    relation = UCSC_RefGene_Group
  )

### GENERATE BARPLOT ###

barplot_relation_genes <- ggplot(
  df_plot_gene,
  aes(
    x = group,
    y = n,
    fill = relation
  )
) +
  geom_bar(
    stat = "identity",
    position = "fill"
  ) +
  ylab("Frequency (%)") +
  xlab("") +
  ggtitle("Relation to Genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 20,
      face = "bold"
    ),
    axis.text.x = element_text(size = 14)
  ) +
  scale_fill_brewer(
    palette = "Dark2"
  )

barplot_relation_genes

ggsave(
  "figures/barplot_relation_genes.png",
  plot = barplot_relation_genes,
  width = 6,
  height = 4,
  dpi = 300
)

### TUMOR CELL CONTENT (TCC) ANALYSIS ###

# Load TCC metadata
tcc <- read.csv(
  "data/EVOFLUx_FFPE_patientData.csv",
  header = TRUE,
  sep = ","
)

# Remove duplicated samples
pca_df$Sample <- rownames(pca_df)

pca_df <- pca_df[
  !pca_df$Sample %in% duplicated_samples,
]

### MATCH TCC METADATA ###

metadata <- tcc[
  grepl("^PBL", tcc$X),
]

pca_df$TCC <- metadata$PURITY_TUMOR_CONSENSUS[
  match(
    pca_df$Sample,
    metadata$X
  )
]

### PCA COLORED BY TCC ###

pca_tcc <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = TCC
  )
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = Sample),
    vjust = -0.5,
    size = 3
  ) +
  scale_color_viridis(
    option = "C",
    direction = 1
  ) +
  labs(
    title = "PCA: PC1 vs PC2 (colored by TCC)",
    x = paste0(
      "PC1 (",
      percent_var[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent_var[2],
      "%)"
    ),
    color = "TCC"
  ) +
  theme_minimal()

pca_tcc

ggsave(
  "figures/pca_tcc.png",
  plot = pca_tcc,
  width = 6,
  height = 4,
  dpi = 300
)

### PCA USING TCC GROUPS ###

pca_df$TCC_group <- ifelse(
  is.na(pca_df$TCC),
  "Missing",
  ifelse(
    pca_df$TCC > 60,
    "High TCC (>60%)",
    "Low TCC (<=60%)"
  )
)

pca_tcc_60 <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = TCC_group
  )
) +
  geom_point(size = 3) +
  geom_text(
    aes(label = Sample),
    vjust = -0.5,
    size = 3
  ) +
  scale_color_manual(
    values = c(
      "High TCC (>60%)" = "#440154FF",
      "Low TCC (<=60%)" = "#21908CFF",
      "Missing" = "gray70"
    )
  ) +
  labs(
    title = "PCA: PC1 vs PC2 (TCC groups)",
    x = paste0(
      "PC1 (",
      percent_var[1],
      "%)"
    ),
    y = paste0(
      "PC2 (",
      percent_var[2],
      "%)"
    ),
    color = "TCC group"
  ) +
  theme_minimal()

pca_tcc_60

ggsave(
  "figures/pca_tcc_60.png",
  plot = pca_tcc_60,
  width = 6,
  height = 4,
  dpi = 300
)

