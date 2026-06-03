
### GO AND KEGG ENRICHMENT ANALYSIS ###

# Required packages
library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg38)
library(pheatmap)

### LOAD EPIC ARRAY ANNOTATION ###

annEPICv1 <- getAnnotation(
  IlluminaHumanMethylationEPICanno.ilm10b4.hg38
)

### DEFINE CpG SETS ###

# Significant differentially methylated CpGs
sig_cpgs <- as.character(sondes.dif)

# Background CpGs
all_cpgs <- as.character(
  rownames(matriu.beta)
)

### GO ENRICHMENT ANALYSIS ###

GO <- gometh(
  sig.cpg = sig_cpgs,
  all.cpg = all_cpgs,
  collection = "GO",
  array.type = "EPIC",
  genomic.features = "ALL"
)

### KEGG ENRICHMENT ANALYSIS ###

kegg <- gometh(
  sig.cpg = sig_cpgs,
  all.cpg = all_cpgs,
  collection = "KEGG",
  array.type = "EPIC",
  genomic.features = "ALL"
)

### FILTER SIGNIFICANT PATHWAYS ###

SigPathways.kegg <- kegg[
  kegg$FDR < 0.05,
]

GO.significant <- GO[
  GO$FDR < 0.05,
]

GO.significant <- GO.significant[
  order(
    GO.significant$FDR,
    decreasing = FALSE
  ),
]

kegg.significant <- kegg[
  kegg$FDR < 0.05,
]

kegg.significant <- kegg.significant[
  order(
    kegg.significant$FDR,
    decreasing = FALSE
  ),
]

### RANK KEGG PATHWAYS BY ENRICHMENT PROPORTION ###

# Calculate DE/N ratio
kegg.significant$prop <- (
  kegg.significant$DE /
    kegg.significant$N
)

# Sort pathways
kegg.rank <- kegg.significant[
  order(
    kegg.significant$prop,
    decreasing = TRUE
  ),
]

kegg.rank

### CLASSIFY CpGs ACCORDING TO DIFFERENTIATION ###

probes_class1 <- rownames(
  umap_df[umap_df$class == 1, ]
)

probes_class2 <- rownames(
  umap_df[umap_df$class == 2, ]
)

### VISUALIZE CLASS 2 CpGs ###

submatrix <- sanes.matrix.ii[
  probes_class2,
]

pheatmap(
  submatrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE
)

### VISUALIZE CLASS 1 CpGs ###

submatrix <- sanes.matrix.ii[
  probes_class1,
]

pheatmap(
  submatrix,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = FALSE
)

### DEFINE HYPERMETHYLATED AND HYPOMETHYLATED PROBES ###

probes_hypermethylated <- probes_class1
probes_hypomethylated <- probes_class2

### GO/KEGG ANALYSIS - HYPERMETHYLATED CpGs ###

GO.hyper <- gometh(
  sig.cpg = probes_hypermethylated,
  all.cpg = all_cpgs,
  collection = "GO",
  array.type = "EPIC",
  genomic.features = "ALL"
)

kegg.hyper <- gometh(
  sig.cpg = probes_hypermethylated,
  all.cpg = all_cpgs,
  collection = "KEGG",
  array.type = "EPIC",
  genomic.features = "ALL"
)

### FILTER SIGNIFICANT HYPERMETHYLATED PATHWAYS ###

GO.hyper.sig <- GO.hyper[
  GO.hyper$FDR < 0.05,
]

GO.hyper.sig <- GO.hyper.sig[
  order(
    GO.hyper.sig$FDR,
    decreasing = FALSE
  ),
]

kegg.hyper.sig <- kegg.hyper[
  kegg.hyper$FDR < 0.05,
]

kegg.hyper.sig <- kegg.hyper.sig[
  order(
    kegg.hyper.sig$FDR,
    decreasing = FALSE
  ),
]

### CALCULATE ENRICHMENT PROPORTION ###

kegg.hyper.sig$prop <- (
  kegg.hyper.sig$DE /
    kegg.hyper.sig$N
)

kegg.hyper.rank <- kegg.hyper.sig[
  order(
    kegg.hyper.sig$prop,
    decreasing = TRUE
  ),
]

### DEFINE HYPOMETHYLATED CpGs ###

probes_hypomethylated <- names(
  deltaBeta[deltaBeta < 0.1]
)

length(probes_hypomethylated)

### GO/KEGG ANALYSIS - HYPOMETHYLATED CpGs ###

GO.hypo <- gometh(
  sig.cpg = probes_hypomethylated,
  all.cpg = all_cpgs,
  collection = "GO",
  array.type = "EPIC",
  genomic.features = "ALL"
)

kegg.hypo <- gometh(
  sig.cpg = probes_hypomethylated,
  all.cpg = all_cpgs,
  collection = "KEGG",
  array.type = "EPIC",
  genomic.features = "ALL"
)

### FILTER SIGNIFICANT HYPOMETHYLATED PATHWAYS ###

GO.hypo.sig <- GO.hypo[
  GO.hypo$FDR < 0.05,
]

GO.hypo.sig <- GO.hypo.sig[
  order(
    GO.hypo.sig$FDR,
    decreasing = FALSE
  ),
]

kegg.hypo.sig <- kegg.hypo[
  kegg.hypo$FDR < 0.05,
]

kegg.hypo.sig <- kegg.hypo.sig[
  order(
    kegg.hypo.sig$FDR,
    decreasing = FALSE
  ),
]

### RANK HYPOMETHYLATED KEGG PATHWAYS ###

kegg.hypo.sig$prop <- (
  kegg.hypo.sig$DE /
    kegg.hypo.sig$N
)

kegg.hypo.rank <- kegg.hypo.sig[
  order(
    kegg.hypo.sig$prop,
    decreasing = TRUE
  ),
]

kegg.hypo.rank

