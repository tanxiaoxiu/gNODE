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
setwd("/lustre/home/acct-clswt/clswt-xiaoxiutan/p3_glv/realdata/data_cdiff/data")
count <- t(read.delim("counts.txt", sep = '\t',header = TRUE, row.names = 1,check.names = FALSE))
biomass <- read.delim("biomass.txt",  sep = '\t', check.names = FALSE)
metadata <- read.delim("metadata.txt",  sep = '\t', check.names = FALSE)

composition <- (count) / rowSums(count)
biomass_mean <-  rowMeans(biomass)
absolute <- composition*biomass_mean

#删除在30%的样本中都不存在的微生物
non_zero_ratio <- colSums(absolute != 0) / nrow(absolute)
filtered_absolute <- absolute[, non_zero_ratio >= 0.3]

data_cdiff <- cbind(metadata,filtered_absolute)
write.table(data_cdiff,file ="data_cdiff.txt",row.names = F,col.names = T, sep = "\t",quote = F)


data_cdiff_filter <- data_cdiff[,c("subjectID","measurementid",
                        "Clostridium-hiranonis","Clostridium-difficile","Proteus-mirabilis",
                        "Clostridium-scindens","Ruminococcus-obeum","Clostridium-ramosum",
                        "Bacteroides-ovatus","Akkermansia-muciniphila","Parabacteroides-distasonis","Bacteroides-fragilis",
                        "Bacteroides-vulgatus","Enterococcus-faecalis","Klebsiella-oxytoca","Lactobacillus-reuteri","Roseburia-hominis","Escherichia-coli")]

colnames(data_cdiff_filter) <- c("Subject","Time",
                   "C.hiranonis","C.difficile","P.mirabilis",
                   "C.scindens","R.obeum","C.ramosum",
                   "B.ovatus","A.muciniphila","P.distasonis","B.fragilis",
                   "B.vulgatus","E.faecalis","K.oxytoca","L.reuteri","R.hominis","E.coli")

data_cdiff_filter[,3:18] <- data_cdiff_filter[,3:18]/(10^8)
microbe_names <- colnames(data_cdiff_filter)[3:18]
write.table(microbe_names, file = "microbe_name.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

data_cdiff_filter$Group <- ifelse(data_cdiff_filter$Time <= 28, "Control", "Case")
data_cdiff_filter <- data_cdiff_filter[, c(names(data_cdiff_filter)[1:2], "Group", names(data_cdiff_filter)[3:ncol(data_cdiff_filter)][names(data_cdiff_filter)[3:ncol(data_cdiff_filter)] != "Group"])]
write.table(data_cdiff_filter,file ="data_cdiff_filter.txt",row.names = F,col.names = T, sep = "\t",quote = F)


#处理成eNODE的输入数据
data_cdiff_filter <- read.delim("/lustre/home/acct-clswt/clswt-xiaoxiutan/p3_glv/realdata/data_cdiff/data/data_cdiff_filter.txt",  sep = '\t', check.names = FALSE)

Con <- data_cdiff_filter %>% filter(Group == "Control")
Con_list_by_subject <- Con %>%
  arrange(Subject, Time) %>%  
  split(.$Subject)  


Con_list <- lapply(Con_list_by_subject, function(df) {
  df <- df[, -c(1:3)]
  mat <- as.matrix(df)
  
  if (nrow(mat) == 13 && ncol(mat) == (16)) {
    rownames(mat) <- 1:13
    return(mat)
  } else {
    return(NULL)  # 如果不符合条件，返回NULL
  }
})
Con_list <- Con_list[!sapply(Con_list, is.null)]
true_initial_state <- simplify2array(Con_list)
saveRDS(Con_list,"Control/Input/true_initial_state.rds")
h5save(file = "Control/Input/true_initial_state.h5", true_initial_state)

##
Case <- data_cdiff_filter %>% filter(Group == "Case")
Case_list_by_subject <- Case %>%
  arrange(Subject, Time) %>%  
  split(.$Subject)  


Case_list <- lapply(Case_list_by_subject, function(df) {
  df <- df[, -c(1:3)]
  mat <- as.matrix(df)
  
  if (nrow(mat) == 13 && ncol(mat) == (16)) {
    rownames(mat) <- 1:13
    return(mat)
  } else {
    return(NULL)  
  }
})
Case_list <- Case_list[!sapply(Case_list, is.null)]
true_initial_state <- simplify2array(Case_list)
saveRDS(Case_list,"Case/Input/true_initial_state.rds")
h5save(file = "Case/Input/true_initial_state.h5", true_initial_state)
