rm(list = ls(all = TRUE)) 
library(MASS)
library(glmnet)
library(Rfast)
library(mgcv)


#write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}


#Data pre-processing
setwd("/lustre/home/acct-clswt/clswt-xiaoxiutan/p3_glv/realdata/data_diet/data_2")
count <- t(read.delim("counts.txt", sep = '\t',header = TRUE, row.names = 1,check.names = FALSE))
biomass <- read.delim("biomass.txt",  sep = '\t', check.names = FALSE)
metadata <- read.delim("metadata.txt",  sep = '\t', check.names = FALSE)

composition <- (count) / rowSums(count)
biomass_mean <-  rowMeans(biomass)
absolute <- composition*biomass_mean

#删除在30%的样本中都不存在的微生物
non_zero_ratio <- colSums(absolute != 0) / nrow(absolute)
filtered_absolute <- absolute[, non_zero_ratio >= 0.3]

data_diet <- cbind(metadata,filtered_absolute)
write.table(data_diet,file ="data_diet.txt",row.names = F,col.names = T, sep = "\t",quote = F)

data_diet_filter <- data_diet[, -c(2, 5:8)]
colnames(data_diet_filter)[1:3] <- c("Sample", "Subject", "Time")
data_diet_filter$Group <- ifelse(data_diet_filter$Time < 30, "High_fiber1", 
                                 ifelse(data_diet_filter$Time > 34 & data_diet_filter$Time <= 49, "Low_fiber", 
                                        ifelse(data_diet_filter$Time > 49, "High_fiber2", NA)))


data_diet_filter <- data_diet_filter[, c(1:3, ncol(data_diet_filter), 4:(ncol(data_diet_filter)-1))]
write.table(data_diet_filter,file ="data_diet_filter.txt",row.names = F,col.names = T, sep = "\t",quote = F)

microbe_names <- colnames(data_diet_filter)[5:17]
write.table(microbe_names, file = "microbe_name.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


#处理成eNODE的输入数据
data_diet_filter <- read.delim("/lustre/home/acct-clswt/clswt-xiaoxiutan/p3_glv/realdata/data_diet/data/data_diet_filter.txt",  sep = '\t', check.names = FALSE)
data_diet_filter <- data_diet_filter[, -1]

#High fiber
High_fiber <- data_diet_filter %>% filter(Group == "High_fiber1")
High_fiber_list_by_subject <- High_fiber %>%
  arrange(Subject, Time) %>%  
  split(.$Subject)  


High_fiber_list <- lapply(High_fiber_list_by_subject, function(df) {
  df <- df[, -c(1:3)]
  mat <- as.matrix(df)
  
  if (nrow(mat) == 25 && ncol(mat) == (13)) {
    rownames(mat) <- 1:25
    return(mat)
  } else {
    return(NULL)  
  }
})
High_fiber_list <- High_fiber_list[!sapply(High_fiber_list, is.null)]
true_initial_state <- simplify2array(High_fiber_list)
saveRDS(High_fiber_list,"High_fiber/Input/true_initial_state.rds")
h5save(file = "High_fiber/Input/true_initial_state.h5", true_initial_state)

##Low fiber
Low_fiber <- data_diet_filter %>% filter(Group == "Low_fiber")
Low_fiber_list_by_subject <- Low_fiber %>%
  arrange(Subject, Time) %>%  
  split(.$Subject)  


Low_fiber_list <- lapply(Low_fiber_list_by_subject, function(df) {
  df <- df[, -c(1:3)]
  mat <- as.matrix(df)
  
  if (nrow(mat) == 15 && ncol(mat) == (13)) {
    rownames(mat) <- 1:15
    return(mat)
  } else {
    return(NULL)  
  }
})
Low_fiber_list <- Low_fiber_list[!sapply(Low_fiber_list, is.null)]
true_initial_state <- simplify2array(Low_fiber_list)
saveRDS(Low_fiber_list,"Low_fiber/Input/true_initial_state.rds")
h5save(file = "Low_fiber/Input/true_initial_state.h5", true_initial_state)

#High fiber2
High_fiber2 <- data_diet_filter %>% filter(Group == "High_fiber2")
High_fiber2_list_by_subject <- High_fiber2 %>%
  arrange(Subject, Time) %>%  
  split(.$Subject)  


High_fiber2_list <- lapply(High_fiber2_list_by_subject, function(df) {
  df <- df[, -c(1:3)]
  mat <- as.matrix(df)
  
  if (nrow(mat) == 14 && ncol(mat) == (13)) {
    rownames(mat) <- 1:14
    return(mat)
  } else {
    return(NULL)  
  }
})
High_fiber2_list <- High_fiber2_list[!sapply(High_fiber2_list, is.null)]
true_initial_state <- simplify2array(High_fiber2_list)
saveRDS(High_fiber2_list,"High_fiber2/Input/true_initial_state.rds")
h5save(file = "High_fiber2/Input/true_initial_state.h5", true_initial_state)
