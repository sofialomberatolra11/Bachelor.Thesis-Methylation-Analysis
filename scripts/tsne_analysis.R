
### t-SNE ANALYSIS ###

# Load required libraries
library(Rtsne)
library(ggplot2)
library(viridis)

# Select only numeric columns
# Exclude Sample and TCC metadata columns
num_cols <- sapply(pca_df, is.numeric)

num_cols["Sample"] <- FALSE
num_cols["TCC"] <- FALSE

# Create dataset for t-SNE
tsne_data <- pca_df[, num_cols]

# Set seed for reproducibility
set.seed(123)

# Perform t-SNE dimensionality reduction
tsne_result <- Rtsne(
  tsne_data,
  dims = 2,
  perplexity = 10
)

# Add t-SNE coordinates to dataframe
pca_df$TSNE1 <- tsne_result$Y[, 1]
pca_df$TSNE2 <- tsne_result$Y[, 2]

# Generate t-SNE plot colored by TCC
tsne_plot <- ggplot(
  pca_df,
  aes(x = TSNE1,
      y = TSNE2,
      color = TCC)
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
    title = "t-SNE Analysis (colored by TCC)",
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "TCC"
  ) +
  theme_minimal()

# Display plot
tsne_plot

# Save figure
ggsave(
  "figures/tsne_tcc.png",
  plot = tsne_plot,
  width = 6,
  height = 4,
  dpi = 300
)
