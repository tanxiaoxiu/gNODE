library(ggplot2)
library(tidyr)
library(tibble)
library(gtable)
library(pheatmap)
library(grid)
library(gridExtra)

setwd("~/gNODE/realdata/data_diet/driver")

summary_data <- read.table(
  file = "summary_degree.txt",
  header = TRUE,        
  sep = "\t",           
  stringsAsFactors = FALSE  
)

l_tb <- summary_data
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Degree")

w_tb <- l_tb %>%
  pivot_wider(id_cols = "Driver",
              names_from = "Group",
              values_from = "Degree",
              values_fill=0) %>%
  mutate(across(-Driver, function(x) x )) %>%
  column_to_rownames("Driver")

p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
g <- p1$gtable


g <- gtable_add_rows(g, heights = unit(0.1, "cm"), pos = 0)
png("Figure10.png", width = 5.5, height = 5.5, units = "in", res = 600)
grid.draw(g)  
dev.off()



