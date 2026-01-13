#Figure5
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(readxl)
library(pheatmap)
library(reshape2)
library(igraph)
library(reticulate)
library(tidyr)

#Figure5A
setwd("~/gNODE/realdata/data_cdiff/trajectory/comparison")
base_path1 <- "~/gNODE/realdata/data_cdiff/trajectory/"

Group <- c("Control","Case")
base_folders <- c("NODE", "pNODE", "gNODE")

path <- paste0(base_path1, "comparison")
setwd(path)
all_rrmse <- readRDS("all_rrmse.rds")
all_rrmse$Group <- factor(all_rrmse$Group, levels = c("Control", "Case"))
all_rrmse$Method <- factor(all_rrmse$Method, levels = c("NODE", "pNODE","gNODE"))

# Calculate summary statistics (mean ± SD)
summary_data <- all_rrmse %>%
  group_by(Group, Method) %>%
  summarise(
    mean_rrmse = mean(rrmse, na.rm = TRUE),
    se_rrmse = sd(rrmse, na.rm = TRUE)/sqrt(n()),
    .groups = 'drop'
  )

# Create the plot
plot2 <- ggplot(summary_data, aes(x = Group, y = mean_rrmse, fill = Method)) +
  geom_col(
    position = position_dodge(width = 0.8),  
    width = 0.7, 
    color = NA,  
    linewidth = 0) +  
  geom_errorbar(
    aes(ymin = mean_rrmse - se_rrmse, ymax = mean_rrmse + se_rrmse),
    position = position_dodge(width = 0.8),  
    width = 0.2,
    color = "black",
    linewidth = 0.8
  ) +
  scale_fill_manual(values = c("#F9C89B","#7CA3B8","#CE8A8D")) +
  xlab("Group") +
  ylab("Relative RMSE") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +  
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    text = element_text(size=28),
    axis.title.x = element_text(size = 30),  
    axis.title.y = element_text(size = 30, margin = margin(r = 15, unit = "pt")), 
    axis.text.x = element_text(size = 28, color = "black"),  
    axis.text.y = element_text(size = 28, color = "black"),
    panel.grid.major = element_line(colour=NA),
    panel.background = element_rect(fill="transparent", colour=NA),
    plot.background = element_rect(fill="transparent", colour=NA),
    panel.grid.minor = element_blank(),
    legend.spacing.y = unit(0.5, 'cm')
  ) +
  guides(fill = guide_legend(byrow = TRUE))

# Save the plot
ggsave(
  filename = "Figure5A.png",
  plot = plot2,
  device = "png",
  dpi = 600,
  units = "in",
  width = 12,
  height = 8
)


#Figure5B
P <- 16
setwd("~/gNODE/realdata/data_cdiff/Control/prediction")

alpha <- read.table("alpha_pre.csv", header = FALSE, sep = "\t")
beta <- read.table("beta_pre.csv", header = FALSE, sep = ",")
microbe_names <- read.table("~/gNODE/realdata/data_cdiff/data/microbe_name.txt", header = FALSE, sep = "\t")

microbe_names <- as.vector(microbe_names$V1)
rownames(beta) <- microbe_names
colnames(beta) <- microbe_names
rownames(alpha) <- microbe_names

p1 <- pheatmap(beta,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(c("#004995", "white", "#85242B"))(100),
               breaks = seq(-1, 1, length.out = 101),  
               center = 0,  
               fontsize_row = 12,        
               fontsize_col = 12,        
               legend = F,
               main = " ",
               show_colnames = TRUE,
               show_rownames = FALSE,
               angle_col = 90)

ggsave(filename="Interaction.png",plot=p1,device="png",dpi=600,units="in",width=3.6,height=4.9)

p2 <- pheatmap(alpha,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(c("#004995", "white", "#85242B"))(100),
               breaks = seq(-1, 1, length.out = 101),  
               center = 0,  
               fontsize_row = 14,        
               fontsize_col = 14,        
               legend = TRUE,
               main = "",
               show_colnames = FALSE)

ggsave(filename="Growth.png",plot=p2,device="png",dpi=600,units="in",width=3,height=4)


#Figure5C
P <- 16
setwd("~/gNODE/realdata/data_cdiff/Case/prediction")
alpha <- read.table("alpha_pre.csv", header = FALSE, sep = "\t")
beta <- read.table("beta_pre.csv", header = FALSE, sep = ",")
microbe_names <- read.table("~/gNODE/realdata/data_cdiff/data/microbe_name.txt", header = FALSE, sep = "\t")
microbe_names <- as.vector(microbe_names$V1)
rownames(beta) <- microbe_names
colnames(beta) <- microbe_names
rownames(alpha) <- microbe_names

p1 <- pheatmap(beta,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(c("#004995", "white", "#85242B"))(100),
               breaks = seq(-1, 1, length.out = 101),  
               center = 0,  
               fontsize_row = 12,        
               fontsize_col = 12,        
               legend = F,
               main = " ",
               show_colnames = TRUE,
               show_rownames = FALSE,
               angle_col = 90)
ggsave(filename="Interaction.png",plot=p1,device="png",dpi=600,units="in",width=3.6,height=4.9)

p2 <- pheatmap(alpha,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               scale = "none",
               color = colorRampPalette(c("#004995", "white", "#85242B"))(100),
               breaks = seq(-1, 1, length.out = 101),  
               center = 0,  
               fontsize_row = 14,        
               fontsize_col = 14,        
               legend = TRUE,
               main = "",
               show_colnames = FALSE)
ggsave(filename="Growth.png",plot=p2,device="png",dpi=600,units="in",width=3,height=4)

#Figure5D
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
