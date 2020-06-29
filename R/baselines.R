# causal inference methods
# 
# implementation details:
# if 'X -> Y' they output cd  = 1
# if 'Y -> X' they output cd = 0
# epsilon: confidence (or score)

library(CAM)
library(Hmisc)

source("R/baselines/Slope/Slope.R")
source("R/baselines/Slope/utilities.R")
source("R/baselines/GRAN_EMD/kernelMethod.R")
source("R/baselines/GRAN_EMD/rkhse.R")
source("R/baselines/Slope/resit/code/startups/startupLINGAM.R", chdir = TRUE)
source("R/baselines/Slope/resit/code/startups/startupICML.R", chdir = TRUE)

SlopeWrap <- function(pair){
  res = SlopeInfo(data.frame(pair), alpha = 1.01)
  cd = NA
  epsilon = NA
  if(res$cd == "->") cd = 1
  else if(res$cd == "<-") cd = 0
  else if(res$cd == "--") cd = NA
  list(epsilon = res$eps, cd = cd)
}

IGCIWrap <- function(pair){
  res = IGCI(data.frame(pair))
  cd = NA
  epsilon = NA
  if(res$cd == "->") cd = 1
  else if(res$cd == "<-") cd = 0
  list(epsilon = res$eps, cd = cd)
  
}


IGCIWrap_G <- function(pair){
  res = IGCI_G(data.frame(pair))
  cd = NA
  epsilon = NA
  if(res$cd == "->") cd = 1
  else if(res$cd == "<-") cd = 0
  list(epsilon = res$eps, cd = cd)
  
}


lingamWrapper = function(t) {
  res = tryCatch({
    r = lingamWrap(t)
  }, error = function(e) {
    print(e)
    return(list(
      B = matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2),
      Adj = matrix(c(F, F, F, F), nrow = 2, ncol = 2)
    ))
  })
  C = 1 * res$Adj
  if (sum(C) != 1) {
    causd = "--"
    p_val = 0
  } else{
    if (C[2, 1] == 1) {
      causd = 0#"<-"
      p_val = res$B[1, 2]
    } else if (C[1, 2] == 1) {
      causd = 1#"->"
      p_val = res$B[2, 1]
    } else{
      causd = "--"
      p_val = 0
    }
  }
  print(res)
  return(list(cd = causd, epsilon = p_val))
  #causd
}

CAMWrapper <- function(X){
 
  res = CAM(X, scoreName = "SEMGAM")
  cd = NA
  if(res$Adj[3] == 1) cd = 1
  else if(res$Adj[2] == 1) cd = 0
  return(list(cd = cd, epsilon = res$Score))
}


ICMLWrapper <- function(X){
  res=ICML(X,
            model = train_gp,
            indtest = indtestHsic,
            output = FALSE)
  cd = NA
  epsilon = NA
  if(res$Cd == "->") cd = 1
  else if(res$Cd == "<-") cd = 0
  return(list(cd = cd, epsilon = res$Eps))
}


# ****synthetic data only****
# excript from original paper: Chen 2014 ; same parameters used in Gaussianity measures 2016
# For EMD we use gaussian kernel with width 1/5*dm, 
# where d_m is the median of distances among all input patterns.

RKHSEWrap <- function(t) {
  res = determine_causal_direction_rkhse(t[,1], t[,2])
  cd = NA
  epsilon = NA
  if(res$relation == "x->y") cd = 1
  else if(res$relation == "y->x") cd = 0
  list(cd = cd, epsilon = res$conf)
}


### test gaussianity measures
###
###

GaussianityWrap <- function(t) {
  res = determine_causal_direction(t[, 1], t[, 2])
  cd = NA
  epsilon = NA
  if(res$relation == "x->y") cd = 1
  else if(res$relation == "y->x") cd = 0
  list(cd = cd, epsilon = res$conf)
}
