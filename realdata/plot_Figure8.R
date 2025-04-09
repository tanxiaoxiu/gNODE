library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(igraph)
library(reticulate)
library(tidyr)


#Case
P <- 16
setwd("~/gNODE/realdata/data_cdiff/Case/prediction")

beta <- read.table("beta_pre.csv", header = FALSE, sep = ",")
microbe_names <- read.table("~/gNODE/realdata/data_cdiff/data/microbe_name.txt", header = FALSE, sep = "\t")

microbe_names <- as.vector(microbe_names$V1)
rownames(beta) <- microbe_names
colnames(beta) <- microbe_names

library(pheatmap)
beta_new <- as.data.frame(matrix(0, nrow = nrow(beta), ncol = ncol(beta), dimnames = list(rownames(beta), colnames(beta)))) 
beta_new["C.difficile",] <- beta["C.difficile",]
beta_new[,"C.difficile"] <- beta[,"C.difficile"]

A <- as.matrix(beta_new)
non_zero_indices <- which(A != 0, arr.ind = TRUE)  
non_zero_elements <- A[non_zero_indices]  
sorted_indices <- order(abs(non_zero_elements), decreasing = TRUE)
top_50_count <- floor(length(sorted_indices) / 2)
top_50_indices <- non_zero_indices[sorted_indices[1:top_50_count], ]
A_filtered <- matrix(0, nrow = nrow(A), ncol = ncol(A))
A_filtered[top_50_indices] <- non_zero_elements[sorted_indices[1:top_50_count]]
rownames(A_filtered) <- rownames(A)
colnames(A_filtered) <- colnames(A)
A_filtered <- t(A_filtered)

g <- graph_from_adjacency_matrix(A_filtered, mode = "directed", weighted = TRUE, diag = FALSE)
E(g)$color <- ifelse(E(g)$weight > 0, "#F9A6A0", "#BCCBE5")
V(g)$label <- rownames(beta)
V(g)$color <- "#FFD5AB"
layout_custom <- layout_in_circle(g)  
center_node <- which(V(g)$name == "C.difficile")  
layout_custom[center_node, ] <- c(0, 0)
other_nodes <- setdiff(1:vcount(g), center_node)
angle <- seq(0, 2*pi, length.out = length(other_nodes) + 1)[-1] 
layout_custom[other_nodes, ] <- cbind(cos(angle), sin(angle))  
width_in_pixels <- 14 * 300
height_in_pixels <- 12 * 300
png(file = "Figure8.png", width = width_in_pixels, height_in_pixels, res = 600)
par(mar=c(0, 0, 0, 0)) 
plot(g, layout = layout_custom,  
     edge.width = abs(E(g)$weight)*10^1.4,
     vertex.label.color = "black", 
     vertex.size = 16,  
     vertex.frame.color = "#FFD5AB",
     edge.arrow.size = 1.5, 
     edge.curved = 0.2, 
     vertex.label.cex = 1.6,
     edge.lty = 1)  
dev.off()



















