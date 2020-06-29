rm(list = ls())

source("R/bqcd.R")
source("R/baselines.R")

n_pairs <- 100
n_size <- 1000

set.seed(3)

path <- "data/ANLSMN_pairs/"
ext_all <- c("AN-s/", "LS-s/", "AN/", "LS/", "MN-U/") #c("LS-s/", "AN-s/", "LS/", "LS-s/", "MN-U/")

do_synthetic_benchmark <- function(method){

  for(ext in ext_all){
    # load ground truth for pairs
    pairs_gt <- as.numeric(read.csv(paste0(path, ext, "pairs_gt.txt"), header = FALSE)$V1)

    method_to_run <- switch(method,
                            "bQCCD" = QCCD,
                            "CAM" = CAMWrapper,
                            "EMD" = RKHSEWrap,
                            "LINGAM" = lingamWrapper,
                            "RESIT" = ICMLWrapper,
                            "GRAN" = GaussianityWrap,
                            "Slope" = SlopeWrap,
                            "IGCI" = IGCIWrap_G,
                            "bQCD_qrf" = QCD_qrf,
                            "bQCD_qnn" = QCD_qnn)

    pair_idx <- rep("", n_pairs)
    eps <- rep(0, n_pairs)
    cds <- rep("--", n_pairs)
    time <- rep(0, n_pairs)
    eps <- rep(0, n_pairs)

    n_corr <- 0
    results <- matrix(0,nrow = n_pairs, ncol = 2)

    for(i in 1:n_pairs) {
      t1 <- Sys.time()
      pair <-  read.table(paste0(path, ext, "pair_",i,".txt"),
                          as.is = TRUE, header = TRUE, sep = ",",
                          row.names = 1)
      plot(pair)
      res <- method_to_run(pair)
      t2 <- Sys.time()
      elapsed <- as.numeric(difftime(t2, t1), units = "secs")
      plot(pair)
      pair_idx[i] <- i
      eps[i] <- res$epsilon
      cds[i] <- res$cd
      time[i] <- elapsed

      # dummy prints for debug
      print(paste0("pair = ", i, ", correct = ", pairs_gt[i] == cds[i],
                   ", accuracy = ", mean(pairs_gt[1:i] == cds[1:i])))
      print("-----")
      print(res)
      results[i,1] <- res$cd

      if (!is.na(res$cd)){
        correct <- (res$cd == pairs_gt[i] )
        n_corr <- n_corr + correct
      }
    }
    ncorrect <- sum(results[1:n_pairs,1]  == pairs_gt[1:n_pairs])
    print(paste("correct =", n_corr, "mean =", ncorrect))

    corr <- rep(0, n_pairs)
    corr[pairs_gt == cds] = 1

    resDF <- data.frame(Correct = corr, Eps = eps, Cds = cds,T = time)
    write.table(resDF, file = paste0("results/",method,"_",substr(ext, 1, nchar(ext) - 1),".tab"),
                row.names = F,
                quote = F,
                sep = "\t")

    acc <- ncorrect/100
    acc
  }
}

methods <- list("bQCD_qnn")#, "CAM", "RESIT", "EMD",  "LINGAM", "Slope")
resultsAll <- lapply(methods, function (m) do_synthetic_benchmark(m))
