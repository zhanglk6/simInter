#if(!requireNamespace("devtools", quietly = TRUE))
#  install.packages("devtools")
#devtools::install_github('zhanglk6/caratINT')
library(caratINT)
source("randomization_methods.R")
source("sample_model.R")
source("test_formula.R")

library(dplyr)
library(sandwich)
library(ranger)
library(np)
library(doSNOW)
library(parallel)
library(xtable)
library(bigsplines)
library(caret)


simulate_model <- function(model_num, N = 800, times = 5000, delta = 0, pi = 0.5, seeds = 439, transformed = FALSE) {
  numCores <- 8
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  
  result <- foreach(i = 1:times, .combine = 'rbind', .packages = c("sandwich","np", "dplyr", "ranger","caratINT")) %dopar% {
    set.seed(seeds + i)
    sample_data <- switch(model_num,
                          {
                            source("randomization_methods.R")
                            source("sample_model.R")
                            source("test_formula.R")
                            sample_model1(N, pi, delta, transformed)
                          },
                          {
                            source("randomization_methods.R")
                            source("sample_model.R")
                            source("test_formula.R")
                            sample_model2(N, pi, delta, transformed)
                          },
                          {
                            source("randomization_methods.R")
                            source("sample_model.R")
                            source("test_formula.R")
                            sample_model3(N, pi, delta, transformed)
                          }
                          ,
                          {
                            source("randomization_methods.R")
                            source("sample_model.R")
                            source("test_formula.R")
                            sample_model4(N, pi, delta, transformed)
                          }
    )
    all_outputs(pi, sample_data)
  }
  stopCluster(cl)
  return(result)
}

simulate_model_small <- function(N = 80, times = 5000, delta = 0, pi = 0.5, seeds = 439) {
  numCores <- 8
  cl <- makeCluster(numCores)
  registerDoSNOW(cl)
  
  result <- foreach(i = 1:times, .combine = 'rbind', .packages = c("sandwich","np", "dplyr", "ranger", "caret","caratINT")) %dopar% {
    set.seed(seeds + i)
    source("randomization_methods.R")
    source("sample_model.R")
    source("test_formula.R")
    sample_data <- sample_model_small(N, pi, delta)
    all_outputs(pi, sample_data)
  }
  stopCluster(cl)
  return(result)  
}

compute_reject_prop <- function(df, threshold = 0.05) {
  return(apply(df, 2, function(x) mean(x < threshold)))
}

generate_table <- function(reject_prop, digits = 2) {
  return(100 * matrix(reject_prop, ncol = 4, byrow = TRUE))
}
digits <- 2

simulate_and_compute <- function(model, N, delta, pi, transformed = FALSE) {
  df <- simulate_model(model, N, times = 5000, delta = delta, pi = pi, transformed = transformed)
  reject_prop <- compute_reject_prop(df)
  table <- generate_table(reject_prop)
  return(table)
}

simulate_and_compute_small <- function(N, delta, pi) {
  df <- simulate_model_small(N, times = 5000, delta = delta, pi = pi)
  reject_prop <- compute_reject_prop(df)
  table <- generate_table(reject_prop)
  return(table)
}

transformed = T
N = 800 # 800
pi = 1 / 2
start_time = Sys.time()

table1_no <- simulate_and_compute(1, N, 0, pi, transformed)
end_time = Sys.time()
end_time - start_time
table1_no

table1_yes <- simulate_and_compute(1, N, 0.5, pi, transformed)
end_time = Sys.time()
end_time - start_time
table1_yes

paper_output <- function(pi = 0.5, N = 200, transformed = FALSE) {
  # Model 1
  table1_no <- simulate_and_compute(1, N, 0, pi, transformed)
  table1_yes <- simulate_and_compute(1, N, 0.5, pi, transformed)
  
  # Model 2
  table2_no <- simulate_and_compute(2, N, 0, pi, transformed)
  table2_yes <- simulate_and_compute(2, N, 0.5, pi, transformed)
  
  # Model 3
  table3_no <- simulate_and_compute(3, N, 0, pi, transformed)
  table3_yes <- simulate_and_compute(3, N, 1.5, pi, transformed)
  
  table4_no <- simulate_and_compute(4, N, 0, pi, transformed)
  table4_yes <- simulate_and_compute(4, N, 1.5, pi, transformed)
  
  table_no <- rbind(table1_no, table2_no, table3_no, table4_no)
  table_yes <- rbind(table1_yes, table2_yes, table3_yes, table4_yes)
  table_together <- cbind(table_no, table_yes)
  
  #xtable(table_together, digits = digits)
  table_together
}

paper_output_small <- function(pi = 0.5) {
  # Model 1
  table1_no <- simulate_and_compute_small(80, 0, pi)
  table1_yes <- simulate_and_compute_small(80, 0.5, pi)
  
  # Model 2
  table2_no <- simulate_and_compute_small(120, 0, pi)
  table2_yes <- simulate_and_compute_small(120, 0.5, pi)
  
  # Model 3
  table3_no <- simulate_and_compute_small(160, 0, pi)
  table3_yes <- simulate_and_compute_small(160, 0.5, pi)
  
  table4_no <- simulate_and_compute_small(200, 0, pi)
  table4_yes <- simulate_and_compute_small(200, 0.5, pi)
  
  table_no <- rbind(table1_no, table2_no, table3_no, table4_no)
  table_yes <- rbind(table1_yes, table2_yes, table3_yes, table4_yes)
  table_together <- cbind(table_no, table_yes)
  
  #xtable(table_together, digits = digits)
  table_together
}

paper_output_trans <- function(pi = 0.5, N = 800) {
  # Model 1
  table1_no <- simulate_and_compute(1, N, 0.5, pi, transformed = F)
  table1_yes <- simulate_and_compute(1, N, 0.5, pi, transformed = T)
  
  # Model 2
  table2_no <- simulate_and_compute(2, N, 0.5, pi, transformed = F)
  table2_yes <- simulate_and_compute(2, N, 0.5, pi, transformed = T)
  
  # Model 3
  table3_no <- simulate_and_compute(3, N, 1.5, pi, transformed = F)
  table3_yes <- simulate_and_compute(3, N, 1.5, pi, transformed = T)
  
  table4_no <- simulate_and_compute(4, N, 1.5, pi, transformed = F)
  table4_yes <- simulate_and_compute(4, N, 1.5, pi, transformed = T)
  
  table_no <- rbind(table1_no, table2_no, table3_no, table4_no)
  table_yes <- rbind(table1_yes, table2_yes, table3_yes, table4_yes)
  table_together <- cbind(table_no, table_yes)
  
  #xtable(table_together, digits = digits)
  table_together
}

####### paper content
tab_1 <- paper_output(1/2, 800)
tab_2 <- paper_output(2/3, 800)

tab_small_1 <- paper_output(1/2, 240)
tab_small_2 <- paper_output_small(1/2)

tab_trans_1 <- paper_output_trans(1/2)
tab_trans_1 %>% round(1)

