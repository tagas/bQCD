# script for generating  benchmark data
# 
rm(list = ls())

source("R/utils/data_generators.R")

# 1. sample pair
# 2. flip a coin to decide the ordering, if coin > 0.5, output X,Y if not Y,X
# 3. save pair in one csv file + ground truth other csv file

n_pairs <- 100 # how many pairs to generate
n_size <- 1000 # sample size of each pair

set.seed(3)

path <- "data/ANLSMN_pairs/"

data_generator <- function(model, n_size){
  path_to_store <- paste0("data/ANLSMN_pairs/", model, "/")
  file.remove(paste0(path_to_store,"pairs_gt.txt"))
  
  func_generator <- switch(model,
                          "AN" = sample_AN,
                          "AN-s" = sample_ANs,
                          "MN-U" = sample_MN_u,
                          "LS" = sample_LS,
                          "LS-s" = sample_LSs
                          )  
  
  for(i in 1:n_pairs){
    
    pair <- func_generator(n_size)
    plot(pair)
    coin <- runif(1, 0, 1)
    
    if(coin < 0.5) {
      pair_out <-  data.frame("x" = pair[,2], "y" = pair[,1])
    } else{
      pair_out <- pair
    }
    pair_gt <- ifelse(coin > 0.5, 1, 0)
    
    write.table(pair_out,  paste0(path_to_store,"pair_",i,".txt"),
                sep = ",", col.names = NA,  qmethod = "double")
    
    write(as.matrix(pair_gt), paste0(path_to_store,"pairs_gt.txt"), append = T)
  }
  
}
# generate pairs for all benchmarks 
models_to_generate <- list("AN", "AN-s", "LS", "LS-s", "MN-U")
lapply(models_to_generate, function (m) data_generator(m, n_size))
