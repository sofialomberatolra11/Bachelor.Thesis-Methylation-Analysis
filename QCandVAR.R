install.packages("data.table")
library(data.table)
data <- fread("/Users/sofialomberatolra/Desktop/TFG IDIBAPS/UnSupAnalisi/Beta.matrix.BIN025.DB8.CM20.LCR90.csv")
head(data)
colnames(data)
data <- as.data.table(data)
length(rownames(data))
cols <- grep("^PBL", colnames(data), value = TRUE)
beta.matrix <- data[, ..cols]
beta.matrix
colnames(beta.matrix)
load("/Users/sofialomberatolra/Desktop/TFG IDIBAPS/UnSupAnalisi/RGSet.samples.RData")
RGSet.all #FILTERED
#Remove first sample

###--HEATMAP--##
#Variance & heatmap
beta.heatmap <- beta.matrix

install.packages("matrixStats")
library(matrixStats)
beta_mat <- as.matrix(beta.heatmap) #matriu numèrica
var.beta.heatmap <- rowVars(beta_mat, na.rm = TRUE) #Calcula la variància per fila de manera molt ràpida

var.beta.heatmap <- apply(beta.heatmap, 1, var, na.rm = TRUE)
taula.variancia <- data.frame(
  Sonda = rownames(beta.matrix),
  Variancia = var.beta.heatmap
)

#1000 first
taula.variancia.ordre <- taula.variancia[order(taula.variancia$Variancia, decreasing = TRUE), ]
var.1000 <- taula.variancia.ordre[1:1000, ]
#Graph
library(pheatmap)
library(RColorBrewer)
valors.beta.heatmap.1000 <- beta.heatmap[rownames(beta.heatmap) %in% var.1000$Sonda, ]
colors <- colorRampPalette(c("blue", "white", "red"))(100)

heatmap.var.beta.1000 <- pheatmap(valors.beta.heatmap.1000,
                                  color = colors,
                                  breaks = seq(0, 1, length.out = 101),  
                                  show_rownames = FALSE,
                                  main = "Heatmap valors beta")
heatmap.var.beta.1000
ggsave("/Users/sofialomberatolra/Desktop/TFG IDIBAPS/UnSupAnalisi/PLOTS/Heatmaps/heatmap.var.beta.1000.png", plot = heatmap.var.beta.1000, width = 8, height = 6, dpi = 120)

#500
var.500 <- taula.variancia.ordre[1:500, ]
#Graph
library(pheatmap)
library(RColorBrewer)
valors.beta.heatmap.500 <- beta.heatmap[rownames(beta.heatmap) %in% var.500$Sonda, ]
colors <- colorRampPalette(c("blue", "white", "red"))(100)

heatmap.var.beta.500 <- pheatmap(valors.beta.heatmap.500,
                                 color = colors,
                                 breaks = seq(0, 1, length.out = 101),  #colors han d'anar de 0 a 1
                                 show_rownames = FALSE,
                                 main = "Heatmap valors beta")
heatmap.var.beta.500
ggsave("/Users/sofialomberatolra/Desktop/TFG IDIBAPS/UnSupAnalisi/PLOTS/Heatmaps/heatmap.var.beta.500.png", plot = heatmap.var.beta.500, width = 8, height = 6, dpi = 120)


#5000
var.5000 <- taula.variancia.ordre[1:5000, ]
#Graph
library(pheatmap)
library(RColorBrewer)
valors.beta.heatmap.5000 <- beta.heatmap[rownames(beta.heatmap) %in% var.5000$Sonda, ]
colors <- colorRampPalette(c("blue", "white", "red"))(100)

heatmap.var.beta.5000 <- pheatmap(valors.beta.heatmap.5000,
                                  color = colors,
                                  breaks = seq(0, 1, length.out = 101),  #colors han d'anar de 0 a 1
                                  show_rownames = FALSE,
                                  main = "Heatmap valors beta")
heatmap.var.beta.5000
ggsave("/Users/sofialomberatolra/Desktop/TFG IDIBAPS/UnSupAnalisi/PLOTS/Heatmaps/heatmap.var.beta.5000.png", plot = heatmap.var.beta.5000, width = 8, height = 6, dpi = 120)