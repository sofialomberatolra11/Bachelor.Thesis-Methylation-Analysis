
### PCA ANALYSIS ###

# Load required libraries
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(carData)
library(car)
library(factoextra)
library(corrplot)
library(PCAtools)
library(viridis)

# Transpose beta matrix to set patients as samples
t_beta <- t(beta.matrix)

# Perform Principal Component Analysis (PCA)
pca_result <- prcomp(t_beta, scale. = TRUE)

# Variance explained by each principal component
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Convert variance to percentage
percent_var <- round(100 * var_explained, 2)

# Display first components
percent_var[1:5]

# Create PCA dataframe
pca_df <- as.data.frame(pca_result$x)
pca_df$Sample <- rownames(pca_df)

### PCA PLOT - UNCOLORED ###

pca_uncolored <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(size = 2) +
  geom_text(aes(label = Sample), vjust = -0.5, size = 3) +
  labs(title = "PCA: PC1 vs PC2",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)")) +
  theme_minimal()

ggsave("figures/pca_uncolored.png",
       plot = pca_uncolored,
       width = 6,
       height = 4,
       dpi = 300)

### LOAD METADATA ###

metadata <- read.csv("data/perSampleQC.csv",
                     header = TRUE,
                     sep = ",")

# Match sample names
pca_df$Sample <- rownames(pca_df)

### PCA COLORED BY BIN SCORE ###

pca_df$BIN.score <- metadata$BIN.score[
  match(pca_df$Sample, metadata$Sample.ID)
]

pca_bins <- ggplot(pca_df,
                   aes(x = PC1,
                       y = PC2,
                       color = BIN.score)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample),
            vjust = -0.5,
            size = 3) +
  scale_color_viridis(option = "C", direction = 1) +
  labs(title = "PCA: PC1 vs PC2 (colored by BIN score)",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"),
       color = "BIN score") +
  theme_minimal()

ggsave("figures/pca_binscore.png",
       plot = pca_bins,
       width = 6,
       height = 4,
       dpi = 300)

### PCA COLORED BY DB SCORE ###

pca_df$DB.score <- metadata$DB.score[
  match(pca_df$Sample, metadata$Sample.ID)
]

pca_DB <- ggplot(pca_df,
                 aes(x = PC1,
                     y = PC2,
                     color = DB.score)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample),
            vjust = -0.5,
            size = 3) +
  scale_color_viridis(option = "D", direction = 1) +
  labs(title = "PCA: PC1 vs PC2 (colored by DB score)",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"),
       color = "DB score") +
  theme_minimal()

ggsave("figures/pca_dbscore.png",
       plot = pca_DB,
       width = 6,
       height = 4,
       dpi = 300)

### PCA COLORED BY CM LOW ###

pca_df$CM.low <- metadata$CM.low[
  match(pca_df$Sample, metadata$Sample.ID)
]

pca_CM.low <- ggplot(pca_df,
                     aes(x = PC1,
                         y = PC2,
                         color = CM.low)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample),
            vjust = -0.5,
            size = 3) +
  scale_color_gradient(low = "blue",
                       high = "red") +
  labs(title = "PCA: PC1 vs PC2 (colored by CM.low)",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"),
       color = "CM low") +
  theme_minimal()

ggsave("figures/pca_cm_low.png",
       plot = pca_CM.low,
       width = 6,
       height = 4,
       dpi = 300)

### PCA COLORED BY CM HIGH ###

pca_df$CM.high <- metadata$CM.high[
  match(pca_df$Sample, metadata$Sample.ID)
]

pca_CM.high <- ggplot(pca_df,
                      aes(x = PC1,
                          y = PC2,
                          color = CM.high)) +
  geom_point(size = 3) +
  geom_text(aes(label = Sample),
            vjust = -0.5,
            size = 3) +
  scale_color_gradient(low = "green",
                       high = "orange") +
  labs(title = "PCA: PC1 vs PC2 (colored by CM.high)",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"),
       color = "CM high") +
  theme_minimal()

ggsave("figures/pca_cm_high.png",
       plot = pca_CM.high,
       width = 6,
       height = 4,
       dpi = 300)

### PCA COLORED BY CM DIFFERENCE ###

pca_df$CM.diff <- metadata$CM.diff[
  match(pca_df$Sample, metadata$Sample.ID)
]

pca_CM.diff <- ggplot(pca_df,
                      aes(x = PC1,
                          y = PC2,
                          color = CM.diff)) +
  geom_point(size = 6) +
  geom_text(aes(label = Sample),
            vjust = -0.5,
            size = 3) +
  scale_color_gradient(low = "blue",
                       high = "yellow") +
  labs(title = "PCA: PC1 vs PC2 (colored by CM difference)",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"),
       color = "CM difference") +
  theme_minimal()

ggsave("figures/pca_cm_difference.png",
       plot = pca_CM.diff,
       width = 6,
       height = 4,
       dpi = 300)

### PCA COLORED BY LCR SCORE ###

pca_df$LCR.score <- metadata$LCR.score[
  match(pca_df$Sample, metadata$Sample.ID)
]

pca_LCR.score <- ggplot(pca_df,
                        aes(x = PC1,
                            y = PC2,
                            color = LCR.score)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample),
            vjust = -0.5,
            size = 3) +
  scale_color_gradient(low = "pink",
                       high = "blue") +
  labs(title = "PCA: PC1 vs PC2 (colored by LCR score)",
       x = paste0("PC1 (", percent_var[1], "%)"),
       y = paste0("PC2 (", percent_var[2], "%)"),
       color = "LCR score") +
  theme_minimal()

ggsave("figures/pca_lcr_score.png",
       plot = pca_LCR.score,
       width = 6,
       height = 4,
       dpi = 300)

