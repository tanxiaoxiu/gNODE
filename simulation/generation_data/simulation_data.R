library(deSolve)
library(tidyverse)

GLV <- function(t, x, parameters) {
  with(as.list(c(x, parameters)), {
    x[x < 10^-5] <- 0
    dxdt <- x * (r + A %*% x)
    list(dxdt)
  })
}


sigma=0.1
sparsity=0 ## 0.2;0.5;0.8

generate_glv_parameters <- function(p, seed) {
  set.seed(seed)
  A <- {
    A_temp <- matrix(rnorm(p * p, mean = 0, sd = sigma), p, p)
    A_temp[sample(1:length(A_temp), round(sparsity * length(A_temp)))] <- 0
    diag(A_temp) <- -1
    A_temp
  }
  
  r <- runif(p)
  
  return(list(A = A, r = r))
}

run_simulation_and_save <- function(s_values, p_values, iter = 10, maxtime = 30, steptime = 1) {
  timepoints <- list(
    t25 = c(0:24),
    t20 = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24),
    t15 = c(0, 1, 2, 3, 4, 6, 7, 9, 11, 13, 15, 17, 19, 21, 24),
    t10 = c(0, 1, 2, 3, 7, 10, 12, 15, 18, 24),
    t5 = c(0, 6, 12, 18, 24),
    t3 = c(0, 12, 24))
  
  
  for (p in p_values) {  
    for (s in s_values) {  
      base_work_dir <- sprintf("~/gNODE/simulation/generation_data/sparsity_0/p%d/s%d", p, s)
      
      for (timepoint_name in names(timepoints)) {  
        timepoint_dir <- file.path(base_work_dir, timepoint_name)
        if (!dir.exists(timepoint_dir)) {
          dir.create(timepoint_dir, recursive = TRUE)
        }
        
        for (iter_idx in 1:iter) { 
          iter_dir <- file.path(timepoint_dir, paste0("iter", iter_idx))
          if (!dir.exists(iter_dir)) {
            dir.create(iter_dir, recursive = TRUE)
          }
          
          parameters <- generate_glv_parameters(p, seed = iter_idx)
          A <- parameters$A
          r <- parameters$r
          
          results <- vector("list", s)
          
          for (subject in 1:s) {  

            set.seed(iter_idx * 100 + subject)
            x0 <- runif(p)
            x0 <- as.matrix(x0)
            
            times <- seq(0, maxtime, by = steptime)
            out <- ode(y = x0, times = times, func = GLV, parms = list(A = A, r = r), method = "ode45")
            
            sparse_indices <- which(out[, "time"] %in% timepoints[[timepoint_name]])
            sparse_out <- out[sparse_indices, , drop = FALSE]
            
            
            new_out_list <- list()
            new_out_sparse_list <- list()
            
            for (n in 1:p) {
              A_new <- A
              r_new <- r
              x0_new <- x0
              x0_new[n] <- 0
              out_new <- ode(y = x0_new, times = times, func = GLV, parms = list(A = A_new, r = r_new), method = "ode45")
            
              sparse_out_new <- out_new[sparse_indices, , drop = FALSE]
              
              new_out_list[[n]] <- out_new
              new_out_sparse_list[[n]] <- sparse_out_new
            }
            
            
            results[[subject]] <- list(
              out = out,
              out_sparse = sparse_out,
              new_out_list = new_out_list,
              new_out_sparse_list = new_out_sparse_list,
              parameter_matrix = list(A = A, r = r)  
            )
          }
          
          saveRDS(results, file = file.path(iter_dir, "results.rds"))
        }
      }
    }
  }
  cat("Simulation complete. Results saved in p/s/t/iter format.\n")
}

#S_values <- c(5, 10, 20, 40)  
#p_values <- c(10)  
#iter <- 10
s_values <- c(5,10,20,40,60) 
#p_values <- c(5,10)
p_values <- c(20,30)
iter <- 10
run_simulation_and_save(s_values = s_values, p_values = p_values, iter = iter)



library(tidyr) 
library(rhdf5)

t_values <- c(5, 10, 15, 20, 25)   
iter_values <- 1:10 

for (p in p_values) {
  for (s in s_values) {
    for (t in t_values) {
      for (iter in iter_values) {
        work_dir <- sprintf("~/gNODE/simulation/generation_data/sparsity_0/p%d/s%d/t%d/iter%d", p, s, t, iter)
        setwd(work_dir)
        
        if (file.exists("results.rds")) {
          results <- readRDS("results.rds")
          
          sparse_out_list <- lapply(results, function(x) x$out_sparse)
          sparse_out_list <- lapply(sparse_out_list, function(x) {
            x <- x[, -which(colnames(x) == "time"), drop = FALSE]
            return(x)
          })
          saveRDS(sparse_out_list, "true_initial_state.rds")
          
          true_initial_state <- simplify2array(sparse_out_list)
          h5save(file = "true_initial_state.h5", true_initial_state)
          
          sparse_new_out_list <- lapply(results, function(simulation) {
            lapply(simulation$new_out_sparse_list, function(data) {
              data[, -which(colnames(data) == "time"), drop = FALSE]
            })
          })
          saveRDS(sparse_new_out_list, "true_new_state.rds")
          
          A <- lapply(results, function(x) x$parameter_matrix$A)
          A <- A[[1]]
          write.table(A, "A_true.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
          
          A_binary <- A
          A_binary[A_binary > 0] <- 1
          write.table(A_binary, "A_true_binary.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
          
          r <- lapply(results, function(x) x$parameter_matrix$r)
          r <- r[[1]]
          write.table(r, "r_true.csv", row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
          
          ##################Score##################

          true_initial_state <- readRDS("true_initial_state.rds")
          true_new_state <- readRDS("true_new_state.rds")
          
          # L2
          calculate_L2_Distances <- function(true_initial_state, true_new_state, s, p) {
            results <- vector("list", s)
            
            for (subject in 1:s) {
              true_initial_state_del <- true_initial_state[[subject]]
              true_new_state_subject <- true_new_state[[subject]]
              
              micro_distance <- numeric(p)
              
              for (n in 1:p) {
                true_new_state_subject_del <- true_new_state_subject[[n]]
                # Calculate L2 distances
                micro_distance[n] <- sqrt(sum((true_initial_state_del[, -n] - true_new_state_subject_del[, -n])^2))
              }
              
              results[[subject]] <- list(micro_distance = micro_distance)
            }
            return(results)
          }
          
          calculate_l2_results <- calculate_L2_Distances(true_initial_state, true_new_state, s, p)
          
          ##################real##################
          # L2
          micro_l2_d_avg <- colMeans(do.call(rbind, lapply(calculate_l2_results, `[[`, "micro_distance")))
          microbes <- paste0("Microbe", 1:p)
          micro_l2_d_avg <- data.frame(Regulate = microbes, Score = micro_l2_d_avg, Group = "Microbe")
          micro_l2_d_avg$Score_normalized <- micro_l2_d_avg$Score / sum(micro_l2_d_avg$Score)
          write.table(micro_l2_d_avg, "score_true_mean_l2.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
          
          ##################subject##################
          micro_l2_d_list <- do.call(rbind, lapply(calculate_l2_results, `[[`, "micro_distance"))
          subjects <- paste0("Subject", 1:s)
          
          micro_l2_normalized <- as.data.frame(sweep(micro_l2_d_list, 1, rowSums(micro_l2_d_list), FUN = "/"))
          colnames(micro_l2_normalized) <- microbes
          micro_l2_normalized <- cbind(Subject = subjects, micro_l2_normalized)
          micro_l2_normalized_long <- pivot_longer(micro_l2_normalized, cols = microbes, names_to = "Regulate", values_to = "Score")
          micro_l2_normalized_long$Group <- "Microbe"
          write.table(micro_l2_normalized_long, "score_normalized_true_subject_l2.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
          
          micro_l2 <- data.frame(micro_l2_d_list)
          colnames(micro_l2) <- microbes
          micro_l2 <- cbind(Subject = subjects, micro_l2)
          micro_l2_long <- pivot_longer(micro_l2, cols = microbes, names_to = "Regulate", values_to = "Score")
          micro_l2_long$Group <- 'Microbe'
          write.table(micro_l2_long,file ="score_true_subject_l2.txt",row.names = F,col.names = TRUE, sep = "\t",quote = F)

          save.image(file = "Save_data.RData")
          
        } else {
          cat("File 'results.rds' does not exist in", work_dir, "\n")
        }
      }
    }
  }
}
