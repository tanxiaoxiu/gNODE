library(ggplot2)
library(tidyr)
library(tibble)
library(gtable)

setwd("~/gNODE/realdata/data_diet/driver")

summary <- read.table(
  file = "summary.txt",
  header = TRUE,        
  sep = "\t",           
  stringsAsFactors = FALSE  
)

l_tb <- summary
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Normalized_Score")

w_tb <- l_tb %>%
  pivot_wider(id_cols = "Driver",
              names_from = "Group",
              values_from = "Normalized_Score",
              values_fill=0) %>%
  mutate(across(-Driver, function(x) x )) %>%
  column_to_rownames("Driver")

p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="Figure9.png",plot=p1,device="png",dpi=600,units="in",width=5,height=6)





