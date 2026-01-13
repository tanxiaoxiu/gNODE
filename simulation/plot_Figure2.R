#Figure2
library(tidyverse)
library(readxl)
library(pheatmap)
library(ggplot2)
library(reshape2)
library(igraph)
library(reticulate)
library(tidyr)
library(cowplot)
library(purrr)
library(data.table)
library(patchwork)

###############setting
base_path <- "~/gNODE/simulation/result1/"
Pi_values <- c(0, 0.2, 0.5, 0.8)
P_values <- c(5,10)
s_values <- c(10,20)
t_values <- c(5,10,15,20,25)
Iter <- c(10)

base_folders <- c("sRidge","pNODE","gNODE")
base_path3 <- paste0(base_path,"comparison_all")

############plot
###鲜亮"#a9d5a3","#9cbedb","#f08f92"
base_path3 <- paste0(base_path,"comparison_all")
for (P in P_values) {
  path <- paste0(base_path3, "/p", P)
  setwd(path)
  all_results <- readRDS("all_pe_r_rmse_mean.rds")
  
  #beta
  beta_r_rmse <- subset(all_results, Group == "beta")
  beta_r_rmse$Method <- factor(beta_r_rmse$Method, levels = c("sRidge","pNODE","gNODE"))
  beta_r_rmse$SubjectN <- factor(beta_r_rmse$SubjectN)
  beta_r_rmse$Timepoint <- as.factor(beta_r_rmse$Timepoint)
  
  pi_labels <- c(
    "0" = "italic(pi)~' ='~0",
    "0.2" = "italic(pi)~' ='~0.2",
    "0.5" = "italic(pi)~' ='~0.5",
    "0.8" = "italic(pi)~' ='~0.8"
  )
  
  subject_labels <- c(
    "10" = "italic(S)~'='~10",
    "20" = "italic(S)~'='~20"
  )
  
  plot_beta_r_rmse <- ggplot(beta_r_rmse, aes(x = Timepoint, y = r_rmse, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(Pi ~ SubjectN,labeller = labeller(Pi = as_labeller(pi_labels,label_parsed),SubjectN = as_labeller(subject_labels,label_parsed))) +
    scale_fill_manual(values = c("#A5D1B0","#7CA3B8","#CE8A8D")) +
    labs(x = expression(italic(T)), y = "Relative RMSE", fill = "Method") +
    theme_minimal() +
    theme(
      strip.background = element_rect(fill = "lightgrey"),
      strip.text = element_text(size = 34,color = "black"),
      axis.title.x = element_text(size = 33,color = "black", margin = margin(t = 19)), 
      axis.title.y = element_text(size = 33,color = "black", margin = margin(r = 18)), 
      axis.text.x = element_text(size = 30,color = "black"),  
      axis.text.y = element_text(size = 30,color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),  
      panel.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
      legend.position = "bottom",  
      legend.title = element_text(size = 33, margin = margin(r = 20)), 
      legend.text = element_text(size = 32, margin = margin(r = 10)),  
      legend.key.size = unit(1, "cm")
    )
  ggsave(filename="Figure2.png",plot=plot_beta_r_rmse,device="png",dpi=600,units="in",width=18,height=20)
  
}


