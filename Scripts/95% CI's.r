#!/bin/R 

#Loading libraries
library(tidyverse)
library(RcppParallel)
library(tidyr)
library(dplyr)
library(multidplyr)

#Defining Function
get_binom_ci <- function(x, n, ci){
  if(n <= 0){
    out <- data.frame(lwr = NA, upr = NA)
  } else {
    out <- binom.test(x, n, ci)$conf.int
    out <- data.frame(lwr = out[1], upr = out[2])
  }
  out
}

#Importing Data
tmp_out <- read.table("confintervals_split.txt", head=F)
tmp_col <- read.table("confintervals_data_names.txt", head=F)
tmp_col[,1] <- as.character(tmp_col[,1])
colnames(tmp_out) <- t(tmp_col)
tmp_names <- as.factor(as.matrix(read.table("Sample_names.txt", head=F)))

cluster <- new_cluster(15)
cluster_copy(cluster, 'get_binom_ci') 

#95% CI
out <- tmp_out %>%
  pivot_longer(cols = -Loc,
               names_to = c('Samp', '.value'),
               names_pattern = "(.*)_T_([_a-z]+)$") %>%
  rowwise %>%
  partition(cluster) %>%
  mutate(CI = get_binom_ci(match, tot_count, ci = 0.95)) %>%
  collect()

save(out, file="confintervals_complete")
