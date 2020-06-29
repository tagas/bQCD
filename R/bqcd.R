# causal inference method
# 
# implementation details:
# if 'X -> Y'  output cd = 1
# if 'Y -> X'  output cd = 0
# epsilon: confidence (or score)
# 
library(quantregForest)
library(rvinecopulib)
library(statmod)
library(qrnn)

## Proper scoring rule for predicted quantiles
quantileScoring <- function(actual, pred, prob = 0.95) {
  mean((as.numeric(actual <= pred) - prob) * (pred - actual))
}

## Quantile Causal discovery method nonparametric copula
QCCD <- function(pair, m=1) {
  
  # to get reproducible jittering results for discreet data
  set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  cop <- bicop(data = u,
               family_set = "tll",
               nonpar_method = "constant")
  
  # deal with discrete data
  pair_scaled <- qnorm(u)
  
  # integrate over quantiles
  if(n < 200){
    uw <- gauss.quad.prob(1)
  } else {
    uw <- gauss.quad.prob(m)
  }
  
  cls <- sapply(uw$nodes, function(uu) {
    u_pred <- cbind(predict(object = cop,
                            newdata = cbind(uu, u[, 2]),
                            what = "hinv2"),
                    predict(object = cop,
                            newdata = cbind(u[, 1], uu),
                            what = "hinv1"))
    
    # marginal and conditional quantiles
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    cond_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], u_pred[, i]))
    
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], cond_q[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
  
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)

  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}


## Quantile Causal discovery method - quantile forest
QCD_qrf <- function(pair, m=3) {
  
  # to get reproducible jittering results for discreet data
  set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  # deal with discrete data
  pair_scaled <- qnorm(u)
  
  # integrate over quantiles
  if(n < 200){
    uw <- gauss.quad.prob(1)
  } else {
    uw <- gauss.quad.prob(m)
  }
  
  colnames(pair_scaled) <- c("x", "y")
  pair_scaled <- as.data.frame(pair_scaled)
  qrf_x_to_y <-  quantregForest(x=as.matrix(pair_scaled[,1]), y=as.matrix(pair_scaled[,2]), nodesize=10,sampsize=50)
  qrf_y_to_x <- quantregForest(x=as.matrix(pair_scaled[,2]), y=as.matrix(pair_scaled[,1]), nodesize=10,sampsize=50)
  
  cls <- sapply(uw$nodes, function(uu) {
    
    pred <- cbind(predict(qrf_y_to_x, as.matrix(pair_scaled[,2]), what=uu),
                  predict(qrf_x_to_y,  as.matrix(pair_scaled[,1]), what=uu))
    
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], pred[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
  
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)
  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}



## Quantile Causal discovery method - Neural net with quantile loss (simoultaneous estimation)
QCD_qnn <- function(pair, m=3) {
  # to get reproducible jittering results for discreet data
  set.seed(0)
  
  n <- nrow(pair)
  # Recover pseudo-observations and estimate the copula non-parametrically
  u <- apply(pair, 2, function(x) rank(x, ties.method = "random")/(n + 1))
  # deal with discrete data
  pair_scaled <- qnorm(u)
  
  # integrate over quantiles
  if(n < 200){
    uw <- gauss.quad.prob(1)
  } else {
    uw <- gauss.quad.prob(m)
  }
  
  if(length(uw$nodes) == 1){
    qnn_x_to_y <- qrnn2.fit(x=as.matrix(pair_scaled[,1]), y=as.matrix(pair_scaled[,2]), tau=uw$nodes,
                            n.hidden=10, n.hidden2=5, n.trials=1,
                            iter.max=1000)
    qnn_y_to_x <- qrnn2.fit(x=as.matrix(pair_scaled[,2]), y=as.matrix(pair_scaled[,1]), tau=uw$nodes,
                            n.hidden=5, n.hidden2=10, n.trials=1,
                            iter.max=1000)
    pred.x_to_y <- qrnn2.predict(as.matrix(pair_scaled[,1]), qnn_x_to_y)
    pred.y_to_x <- qrnn2.predict(as.matrix(pair_scaled[,2]), qnn_y_to_x)
  }
  else{
    qnn_x_to_y <- mcqrnn.fit(x=as.matrix(pair_scaled[,1]), y=as.matrix(pair_scaled[,2]), tau=uw$nodes,
                             n.hidden=10, n.hidden2=5, n.trials=1,
                             iter.max=1000)
    qnn_y_to_x <- mcqrnn.fit(x=as.matrix(pair_scaled[,2]), y=as.matrix(pair_scaled[,1]), tau=uw$nodes,
                             n.hidden=5, n.hidden2=10, n.trials=1,
                             iter.max=1000)
    pred.x_to_y <- mcqrnn.predict(as.matrix(pair_scaled[,1]), qnn_x_to_y)
    pred.y_to_x <- mcqrnn.predict(as.matrix(pair_scaled[,2]), qnn_y_to_x)
  }
  
  cls <- sapply(uw$nodes, function(uu) {
    
    pred <- cbind(pred.y_to_x[ ,which(uw$nodes==uu)],
                  pred.x_to_y[ ,which(uw$nodes==uu)])
    
    marg_q <- sapply(1:2, function(i) quantile(pair_scaled[,i], uu))
    # code lengths
    cl_marginal <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], marg_q[i], uu))
    cl_conditional <- sapply(1:2, function(i)
      quantileScoring(pair_scaled[, i], pred[, i], uu))
    
    c(cl_marginal, cl_conditional)
  })
  
  sel <- !apply(is.na(cls), 2, any)
  uw$weights <- uw$weights[sel] / sum(uw$weights[sel])
  cls <- apply(na.omit(t(cls)) * uw$weights, 2, sum)
  
  dx_to_y <- (cls[1] + cls[4])/sum(cls[1:2])
  dy_to_x <- (cls[2] + cls[3])/sum(cls[1:2])
  
  cd <- ifelse(dy_to_x > dx_to_y, 1, 0)
  epsilon <-  (-dx_to_y + dy_to_x )
  
  return(list(cd = cd, epsilon = epsilon))
}


# wrapper used for the real data pairs
QCD_wrap <- function(X, Y, type="QCCD", m=7){
  method_to_run <- switch(type,
                          "QCCD" = QCCD,
                          "QCD_qnn" = QCD_qnn,
                          "QCD_qrf" = QCD_qrf)
  
  res = method_to_run(cbind(X,Y), m)
  if(!is.na(res$cd)) {
    cd = ifelse(res$cd == 1, "->", "<-")
  } else{
    cd = "--"
  }
  list(cd = cd, eps = res$eps)
}