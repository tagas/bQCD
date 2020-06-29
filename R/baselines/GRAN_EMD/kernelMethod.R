library(energy)
source("R/baselines/Npres_Fucntions.R")

##
# These functions compute the kernel matrices needed by the method to determine the causal direction
#
# Author: Daniel Hernández Lobato, Jun 2014
#

# The matrix x contains in each row one observation.
# Some kenrel functions have parameters such as the Gaussian kernel.
# Tau is a regularization parameter that has to be fixed somehow.

##
# The following functions compute the kernel matrices needed for the method
#

gaussian_kernel <- function(x, sigma) {
  exp(-sigma * as.matrix(dist(x))^2)
}

linear_kernel <- function(x) {
  x %*% t(x)
}

gaussian_kernel_XY <- function(x, y, sigma) {
  squared_distances <- (matrix(diag(x %*% t(x)), nrow(x),nrow(y)) - x %*% t(y) - 
                          t(y %*% t(x)) +  matrix(diag(y %*% t(y)), nrow(x), nrow(y), byrow = TRUE))
  exp(-sigma * squared_distances)
}

linear_kernel_XY <- function(x,y) {
  x %*% t(y)
}

##
# The following functions compute the kernel matrix of the residuals with and without centering the data
#

gaussian_kernel_residuals <- function(x, y, sigma, tau) {
  
  Kxx <- gaussian_kernel(x, sigma)
  Kyy <- gaussian_kernel(y, sigma)
  Kxy <- gaussian_kernel_XY(x, y, sigma)
  Kyx <- t(Kxy)
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxxV <- Kxx %*% V
  KxxVKyy <- KxxV %*% Kyy
  KxxVKyyVKxx <- KxxVKyy %*% t(KxxV)
  
  Kyy - t(KxxVKyy) - KxxVKyy + KxxVKyyVKxx
}

linear_kernel_residuals <- function(x, y, tau) {
  
  Kxx <- linear_kernel(x)
  Kyy <- linear_kernel(y)
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxxV <- Kxx %*% V
  KxxVKyy <- KxxV %*% Kyy
  KxxVKyyVKxx <- KxxVKyy %*% t(KxxV)
  
  Kyy - t(KxxVKyy) - KxxVKyy + KxxVKyyVKxx
}

linear_kernel_residuals_centered_data <- function(x, y, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- linear_kernel(x)
  Kyy <- linear_kernel(y)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxxV <- Kxx %*% V
  KxxVKyy <- KxxV %*% Kyy
  KxxVKyyVKxx <- KxxVKyy %*% t(KxxV)
  
  Kyy - t(KxxVKyy) - KxxVKyy + KxxVKyyVKxx
}

gaussian_kernel_residuals_centered_data <- function(x, y, sigma, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- gaussian_kernel(x, sigma)
  Kyy <- gaussian_kernel(y, sigma)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxxV <- Kxx %*% V
  KxxVKyy <- KxxV %*% Kyy
  KxxVKyyVKxx <- KxxVKyy %*% t(KxxV)
  
  Kyy - t(KxxVKyy) - KxxVKyy + KxxVKyyVKxx
}

linear_kernel_centered_residuals_centered_data <- function(x, y, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- linear_kernel(x)
  Kyy <- linear_kernel(y)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxxV <- Kxx %*% V
  KxxVKyy <- KxxV %*% Kyy
  KxxVKyyVKxx <- KxxVKyy %*% t(KxxV)
  
  Kee <- Kyy - t(KxxVKyy) - KxxVKyy + KxxVKyyVKxx
  Kee - ones %*% Kee - Kee %*% ones + ones %*% Kee %*% ones
}

gaussian_kernel_centered_residuals_centered_data <- function(x, y, sigma, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- gaussian_kernel(x, sigma)
  Kyy <- gaussian_kernel(y, sigma)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxxV <- Kxx %*% V
  KxxVKyy <- KxxV %*% Kyy
  KxxVKyyVKxx <- KxxVKyy %*% t(KxxV)
  
  Kee <- Kyy - t(KxxVKyy) - KxxVKyy + KxxVKyyVKxx
  Kee - ones %*% Kee - Kee %*% ones + ones %*% Kee %*% ones
}

##
# The following functions compute the whitened residuals with a centering of the data
#

linear_kernel_whitened_residuals_centered_data <- function(x, y, tau) {
  
  Kee <- linear_kernel_centered_residuals_centered_data(x, y, tau)
  ret <- eigen(Kee, symmetric = TRUE)
  
  list(residuals = sqrt(nrow(x)) * Re(ret$vectors), values = abs(Re(ret$values) / nrow(x)))
}

gaussian_kernel_whitened_residuals_centered_data <- function(x, y, sigma, tau) {
  
  Kee <- gaussian_kernel_centered_residuals_centered_data(x, y, sigma, tau)
  ret <- eigen(Kee, symmetric = TRUE)
  
  list(residuals = sqrt(nrow(x)) * Re(ret$vectors), values = abs(Re(ret$values) / nrow(x)))
}

##
# The following functions compute the squared prediction error on new data
#

evaluate_squared_error_new_data_linear_kernel_centering <- function(x, y, xnew, ynew, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- linear_kernel(x)
  Kyy <- linear_kernel(y)
  
  Kynewynew <- linear_kernel(ynew)
  Kynewy <- linear_kernel_XY(ynew, y)
  Kyynew <- t(Kynewy)
  Kxnewx <- linear_kernel_XY(xnew, x)
  
  M <- matrix(1 / nrow(x), nrow(xnew), nrow(x))
  
  Kynewynew <- Kynewynew - Kynewy %*% t(M) - M %*% Kyynew + M %*% Kyy %*% t(M)
  Kxnewx <- Kxnewx - M %*% Kxx - Kxnewx %*% ones + M %*% Kxx %*% ones
  Kynewy <- Kynewy - M %*% Kyy - Kynewy %*% ones + M %*% Kyy %*% ones
  Kyynew <- t(Kynewy)
  Kxxnew <- t(Kxnewx)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxnewxV <- Kxnewx %*% V
  KxnewxVKyy <- KxnewxV %*% Kyy
  KxnewxVKyyVKxxnew <- KxnewxVKyy %*% t(KxnewxV)
  KxnewxVKyynew <- KxnewxV %*% Kyynew
  
  E <- Kynewynew - KxnewxVKyynew - t(KxnewxVKyynew) + KxnewxVKyyVKxxnew
  sum(diag(E))
}

evaluate_squared_error_new_data_gaussian_kernel_centering <- function(x, y, xnew, ynew, sigma, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- gaussian_kernel(x, sigma)
  Kyy <- gaussian_kernel(y, sigma)
  
  Kynewynew <- gaussian_kernel(ynew, sigma)
  Kynewy <- gaussian_kernel_XY(ynew, y, sigma)
  Kyynew <- t(Kynewy)
  Kxnewx <- gaussian_kernel_XY(xnew, x, sigma)
  
  M <- matrix(1 / nrow(x), nrow(xnew), nrow(x))
  
  Kynewynew <- Kynewynew - Kynewy %*% t(M) - M %*% Kyynew + M %*% Kyy %*% t(M)
  Kxnewx <- Kxnewx - M %*% Kxx - Kxnewx %*% ones + M %*% Kxx %*% ones
  Kynewy <- Kynewy - M %*% Kyy - Kynewy %*% ones + M %*% Kyy %*% ones
  Kyynew <- t(Kynewy)
  Kxxnew <- t(Kxnewx)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxnewxV <- Kxnewx %*% V
  KxnewxVKyy <- KxnewxV %*% Kyy
  KxnewxVKyyVKxxnew <- KxnewxVKyy %*% t(KxnewxV)
  KxnewxVKyynew <- KxnewxV %*% Kyynew
  
  E <- Kynewynew - KxnewxVKyynew - t(KxnewxVKyynew) + KxnewxVKyyVKxxnew
  sum(diag(E))
}

##
# These functions evaluate the explained variance
#

evaluate_explained_variance_new_data_linear_kernel_centering <- function(x, y, xnew, ynew, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- linear_kernel(x)
  Kyy <- linear_kernel(y)
  
  Kynewynew <- linear_kernel(ynew)
  Kynewy <- linear_kernel_XY(ynew, y)
  Kyynew <- t(Kynewy)
  Kxnewx <- linear_kernel_XY(xnew, x)
  
  M <- matrix(1 / nrow(x), nrow(xnew), nrow(x))
  
  Kynewynew <- Kynewynew - Kynewy %*% t(M) - M %*% Kyynew + M %*% Kyy %*% t(M)
  Kxnewx <- Kxnewx - M %*% Kxx - Kxnewx %*% ones + M %*% Kxx %*% ones
  Kynewy <- Kynewy - M %*% Kyy - Kynewy %*% ones + M %*% Kyy %*% ones
  Kyynew <- t(Kynewy)
  Kxxnew <- t(Kxnewx)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxnewxV <- Kxnewx %*% V
  KxnewxVKyy <- KxnewxV %*% Kyy
  KxnewxVKyyVKxxnew <- KxnewxVKyy %*% t(KxnewxV)
  KxnewxVKyynew <- KxnewxV %*% Kyynew
  
  E <- Kynewynew - KxnewxVKyynew - t(KxnewxVKyynew) + KxnewxVKyyVKxxnew
  1 - sum(diag(E)) / sum(diag(Kynewynew))
  
}

evaluate_explained_variance_new_data_gaussian_kernel_centering <- function(x, y, xnew, ynew, sigma, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- gaussian_kernel(x, sigma)
  Kyy <- gaussian_kernel(y, sigma)
  
  Kynewynew <- gaussian_kernel(ynew, sigma)
  Kynewy <- gaussian_kernel_XY(ynew, y, sigma)
  Kyynew <- t(Kynewy)
  Kxnewx <- gaussian_kernel_XY(xnew, x, sigma)
  
  M <- matrix(1 / nrow(x), nrow(xnew), nrow(x))
  
  Kynewynew <- Kynewynew - Kynewy %*% t(M) - M %*% Kyynew + M %*% Kyy %*% t(M)
  Kxnewx <- Kxnewx - M %*% Kxx - Kxnewx %*% ones + M %*% Kxx %*% ones
  Kynewy <- Kynewy - M %*% Kyy - Kynewy %*% ones + M %*% Kyy %*% ones
  Kyynew <- t(Kynewy)
  Kxxnew <- t(Kxnewx)
  
  Kxx <- Kxx - ones %*% Kxx - Kxx %*% ones + ones %*% Kxx %*% ones
  Kyy <- Kyy - ones %*% Kyy - Kyy %*% ones + ones %*% Kyy %*% ones
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxnewxV <- Kxnewx %*% V
  KxnewxVKyy <- KxnewxV %*% Kyy
  KxnewxVKyyVKxxnew <- KxnewxVKyy %*% t(KxnewxV)
  KxnewxVKyynew <- KxnewxV %*% Kyynew
  
  E <- Kynewynew - KxnewxVKyynew - t(KxnewxVKyynew) + KxnewxVKyyVKxxnew
  1 - sum(diag(E)) / sum(diag(Kynewynew))
}

evaluate_explained_variance_new_data_linear_kernel_centering_fast <- function(x, y, xnew, ynew, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- linear_kernel(x)
  Kyy <- linear_kernel(y)
  
  Kynewynew <- linear_kernel(ynew)
  Kynewy <- linear_kernel_XY(ynew, y)
  Kyynew <- t(Kynewy)
  Kxnewx <- linear_kernel_XY(xnew, x)
  
  M <- matrix(1 / nrow(x), nrow(xnew), nrow(x))
  
  M1 <- matrix(Kynewy %*% M[ 1, ], nrow(Kynewynew), ncol(Kynewynew))
  Kynewynew <- Kynewynew - M1 - t(M1) + sum(Kyy) / nrow(x)^2
  
  M1 <- matrix(M[ 1, ] %*% Kxx, nrow(Kxnewx), ncol(Kxnewx), byrow = TRUE)
  M2 <- matrix(Kxnewx %*% ones[ 1, ], nrow(Kxnewx), ncol(Kxnewx))
  Kxnewx <- Kxnewx - M1 - M2 + sum(Kxx) / nrow(x)^2
  
  M1 <- matrix(M[ 1, ] %*% Kyy, nrow(Kynewy), ncol(Kynewy), byrow = TRUE)
  M2 <- matrix(Kynewy %*% ones[ 1, ], nrow(Kynewy), ncol(Kynewy))
  Kynewy <- Kynewy - M1 - M2 + sum(Kyy) / nrow(x)^2
  
  Kyynew <- t(Kynewy)
  Kxxnew <- t(Kxnewx)
  
  M1 <- matrix(ones[ 1, ] %*% Kxx, nrow(Kxx), ncol(Kxx), byrow = TRUE)
  Kxx <- Kxx - M1 - t(M1) + sum(Kxx) / nrow(x)^2
  
  M1 <- matrix(ones[ 1, ] %*% Kyy, nrow(Kyy), ncol(Kyy), byrow = TRUE)
  Kyy <- Kyy - M1 - t(M1) + sum(Kyy) / nrow(x)^2
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxnewxV <- Kxnewx %*% V
  KxnewxVKyy <- KxnewxV %*% Kyy
  KxnewxVKyyVKxxnew <- KxnewxVKyy %*% t(KxnewxV)
  KxnewxVKyynew <- KxnewxV %*% Kyynew
  
  E <- Kynewynew - KxnewxVKyynew - t(KxnewxVKyynew) + KxnewxVKyyVKxxnew
  1 - sum(diag(E)) / sum(diag(Kynewynew))
}


evaluate_explained_variance_new_data_gaussian_kernel_centering_fast <- function(x, y, xnew, ynew, sigma, tau) {
  
  ones <- matrix(1 / nrow(x), nrow(x), nrow(x))
  
  Kxx <- gaussian_kernel(x, sigma)
  Kyy <- gaussian_kernel(y, sigma)
  
  Kynewynew <- gaussian_kernel(ynew, sigma)
  Kynewy <- gaussian_kernel_XY(ynew, y, sigma)
  Kyynew <- t(Kynewy)
  Kxnewx <- gaussian_kernel_XY(xnew, x, sigma)
  
  M <- matrix(1 / nrow(x), nrow(xnew), nrow(x))
  
  M1 <- matrix(Kynewy %*% M[ 1, ], nrow(Kynewynew), ncol(Kynewynew))
  Kynewynew <- Kynewynew - M1 - t(M1) + sum(Kyy) / nrow(x)^2
  
  M1 <- matrix(M[ 1, ] %*% Kxx, nrow(Kxnewx), ncol(Kxnewx), byrow = TRUE)
  M2 <- matrix(Kxnewx %*% ones[ 1, ], nrow(Kxnewx), ncol(Kxnewx))
  Kxnewx <- Kxnewx - M1 - M2 + sum(Kxx) / nrow(x)^2
  
  M1 <- matrix(M[ 1, ] %*% Kyy, nrow(Kynewy), ncol(Kynewy), byrow = TRUE)
  M2 <- matrix(Kynewy %*% ones[ 1, ], nrow(Kynewy), ncol(Kynewy))
  Kynewy <- Kynewy - M1 - M2 + sum(Kyy) / nrow(x)^2
  
  Kyynew <- t(Kynewy)
  Kxxnew <- t(Kxnewx)
  
  M1 <- matrix(ones[ 1, ] %*% Kxx, nrow(Kxx), ncol(Kxx), byrow = TRUE)
  Kxx <- Kxx - M1 - t(M1) + sum(Kxx) / nrow(x)^2
  
  M1 <- matrix(ones[ 1, ] %*% Kyy, nrow(Kyy), ncol(Kyy), byrow = TRUE)
  Kyy <- Kyy - M1 - t(M1) + sum(Kyy) / nrow(x)^2
  
  V <- chol2inv(chol(Kxx + diag(nrow(Kxx)) * tau))
  
  KxnewxV <- Kxnewx %*% V
  KxnewxVKyy <- KxnewxV %*% Kyy
  KxnewxVKyyVKxxnew <- KxnewxVKyy %*% t(KxnewxV)
  KxnewxVKyynew <- KxnewxV %*% Kyynew
  
  E <- Kynewynew - KxnewxVKyynew - t(KxnewxVKyynew) + KxnewxVKyyVKxxnew
  1 - sum(diag(E)) / sum(diag(Kynewynew))
}

##
# The following functions tune the hyper-parameters using 10-fold-cv. It uses the percenage of explained variance as 
# a metric and a grid search
#

tune_parameters_gaussian_kernel <- function(x, y, nfolds = 10, sigma = exp(seq(-10, 3, length = 10)), tau = exp(seq(-7, 7, length = 10))) {
  
  e <- matrix(0, length(sigma), length(tau))
  
  split <- sample(cut(1 : nrow(x), nfolds))
  
  for (i in 1 : length(sigma)) {
    for (j in 1 : length(tau)) {
      for (k in 1 : length(levels(split))) {
        
        # We randomly partition the data into train and validation
        
        xtrain <- x[ split != levels(split)[ k ],, drop = FALSE ]
        ytrain <- y[ split != levels(split)[ k ],, drop = FALSE ]
        xval <- x[ split == levels(split)[ k ],, drop = FALSE ]
        yval <- y[ split == levels(split)[ k ],, drop = FALSE ]
        
        # We evaluate the explained variance
        
        e[ i, j ] <- e[ i, j ] + evaluate_explained_variance_new_data_gaussian_kernel_centering_fast(xtrain, ytrain, 
                                                                                                     xval, yval, sigma[ i ], tau[ j ]) 
      }
      
      e[ i, j ] <- e[ i, j ] / k
    }
    cat(".")
  }
  
  cat("\n")
  
  best <- -Inf
  bestSigma <- 0
  bestTau <- 0
  
  for (i in 1 : length(sigma))
    for (j in 1 : length(tau)) {
      if (best < e[ i, j]) {
        best <- e[ i ,j ]
        bestSigma <- sigma[ i ]
        bestTau <- tau[ j ]
      }
    }
  
  list(expVar = best, sigma = bestSigma, tau = bestTau)
}

tune_parameters_linear_kernel <- function(x, y, nfolds = 10, tau = exp(seq(-7, 7, length = 10))) {
  
  e <- rep(0, length(tau))
  
  split <- sample(cut(1 : nrow(x), nfolds))
  
  for (j in 1 : length(tau)) {
    for (k in 1 : length(levels(split))) {
      
      # We randomly partition the data into train and validation
      
      xtrain <- x[ split != levels(split)[ k ],, drop = FALSE ]
      ytrain <- y[ split != levels(split)[ k ],, drop = FALSE ]
      xval <- x[ split == levels(split)[ k ],, drop = FALSE ]
      yval <- y[ split == levels(split)[ k ],, drop = FALSE ]
      
      # We evaluate the explained variance
      
      e[ j ] <- e[ j ] + evaluate_explained_variance_new_data_linear_kernel_centering_fast(xtrain, ytrain, xval, yval, tau[ j ]) 
    }
    
    e[ j ] <- e[ j ] / k
  }
  
  best <- -Inf
  bestTau <- 0
  
  for (j in 1 : length(tau)) {
    if (best < e[ j]) {
      best <- e[ j ]
      bestTau <- tau[ j ]
    }
  }
  
  list(expVar = best, tau = bestTau)
}

tune_parameters_gaussian_kernel_heuristic <- function(x, y, nfolds = 10, tau = exp(seq(-7, 7, length = 10))) {
  
  e <- rep(0, length(tau))
  
  split <- sample(cut(1 : nrow(x), nfolds))
  
  distances <- as.matrix(dist(x))
  distances[ lower.tri(distances) ] <- 0
  sigma <- 1 / (2 * median(distances[ distances > 0 ])) 
  
  for (j in 1 : length(tau)) {
    for (k in 1 : length(levels(split))) {
      
      # We randomly partition the data into train and validation
      
      xtrain <- x[ split != levels(split)[ k ],, drop = FALSE ]
      ytrain <- y[ split != levels(split)[ k ],, drop = FALSE ]
      xval <- x[ split == levels(split)[ k ],, drop = FALSE ]
      yval <- y[ split == levels(split)[ k ],, drop = FALSE ]
      
      # We evaluate the explained variance
      
      e[ j ] <- e[ j ] + evaluate_explained_variance_new_data_gaussian_kernel_centering_fast(xtrain, ytrain, 
                                                                                             xval, yval, sigma, tau[ j ]) 
    }
    
    e[ j ] <- e[ j ] / k
  }
  
  best <- -Inf
  bestTau <- 0
  
  for (j in 1 : length(tau)) 
    if (best < e[ j ]) {
      best <- e[ j ]
      bestTau <- tau[ j ]
    }
  
  list(expVar = best, sigma = sigma, tau = bestTau)
}

tune_parameters_gaussian_kernel_heuristic_fixed_tau <- function(x, y) {
  
  distances_x <- as.matrix(dist(x))
  distances_x[ lower.tri(distances_x) ] <- 0
  distances_y <- as.matrix(dist(y))
  distances_y[ lower.tri(distances_y) ] <- 0
  
  sigma <- 1 / (2 * median(c(distances_x[ distances_x > 0 ], distances_y[ distances_y > 0 ]))) 
  tau <- 1e-3
  
  list(sigma = sigma, tau = 1e-3)
}


##
# This function is used to make that both x and y have the same distribution. For this x is transformed to have the same
# distribution as y
#

transformDataX <- function(x, y) {
  
  fvalues <- ecdf(y)(y)
  xvalues <- y
  quantiles_y <- approxfun(fvalues, xvalues)
  cumulative_x <- ecdf(x)
  fvalues <- ecdf(x)(x)
  xvalues <- x
  quantiles_x <- approxfun(fvalues, xvalues)
  cumulative_y <- ecdf(y)
  
  transform <- function(x) quantiles_y(cumulative_x(x))
  itransform <- function(x) quantiles_x(cumulative_y(x))
  x <- transform(x)
  
  # This avoids NAs due to repeated elements
  
  x[ is.na(x) ] <- min(y)
  
  list(x = matrix(x, length(x), 1), y = y, transform = transform, itransform = itransform)
}

##
# We compute the some required values to make a decision
#

method_computations <- function(x, y) {
  
  # We transform x to have the same distribution as y's
  
  x_orig <- x
  y_orig <- y
  
  x <- transformDataX(x, y)$x
  
  ret <- tune_parameters_gaussian_kernel(x, y)
  
  sigma <- ret$sigma
  tau <- ret$tau
  
  resultForward <- gaussian_kernel_whitened_residuals_centered_data(x, y, sigma, tau)
  
  ret <- tune_parameters_gaussian_kernel(y, x)
  
  sigma <- ret$sigma
  tau <- ret$tau
  
  resultBackward <- gaussian_kernel_whitened_residuals_centered_data(y, x, sigma, tau)
  
  # We do the same but chaning y and x (These results will end with the letter R)
  
  y <- x_orig
  x <- y_orig
  
  x <- transformDataX(x, y)$x
  
  ret <- tune_parameters_gaussian_kernel(x, y)
  
  sigma <- ret$sigma
  tau <- ret$tau
  
  resultForwardR <- gaussian_kernel_whitened_residuals_centered_data(x, y, sigma, tau)
  
  ret <- tune_parameters_gaussian_kernel(y, x)
  
  sigma <- ret$sigma
  tau <- ret$tau
  
  resultBackwardR <- gaussian_kernel_whitened_residuals_centered_data(y, x, sigma, tau)
  
  toSelectForward <- which(cumsum(resultForward$values) / sum(resultForward$values) > 0.9)[ 1 ]
  toSelectBackward <- which(cumsum(resultBackward$values) / sum(resultBackward$values) > 0.9)[ 1 ]
  toSelectForwardR <- which(cumsum(resultForwardR$values) / sum(resultForwardR$values) > 0.9)[ 1 ]
  toSelectBackwardR <- which(cumsum(resultBackwardR$values) / sum(resultBackwardR$values) > 0.9)[ 1 ]
  
  toSelectForward <- toSelectBackward <- toSelectBackwardR <- toSelectForwardR <- 1
  
  residualsForward <- resultForward$residuals[ , 1 : toSelectForward, drop = FALSE ] 
  residualsBackward <- resultBackward$residuals[ , 1 : toSelectBackward, drop = FALSE ] 
  residualsForwardR <- resultForwardR$residuals[ , 1 : toSelectForwardR, drop = FALSE ] 
  residualsBackwardR <- resultBackwardR$residuals[ , 1 : toSelectBackwardR, drop = FALSE ] 
  
  # We fit a Gamma distribution to the emprical distribution of the statistic to approximate the p-values
  # We use 1e5 samples from the distribution of the statistic under the null hypothesis
  
  if (toSelectForward  == 1) 
    testResultForward <- std.mvnorm.e(c(residualsForward))
  else 
    testResultForward <-std.mvnorm.e(residualsForward)
  
  if (toSelectBackward == 1) 
    testResultBackward <- std.mvnorm.e(c(residualsBackward))
  else
    testResultBackward <- std.mvnorm.e(residualsBackward)
  
  if (toSelectForwardR == 1) 
    testResultForwardR <- std.mvnorm.e(c(residualsForwardR))
  else
    testResultForwardR <- std.mvnorm.e(residualsForwardR)
  
  if (toSelectBackwardR == 1) 
    testResultBackwardR <- std.mvnorm.e(c(residualsBackwardR))
  else
    testResultBackwardR <- std.mvnorm.e(residualsBackwardR)
  
  # We return all computations, including approximate p-values
  
  list(testResultForward = testResultForward, 
       testResultBackward = testResultBackward,
       testResultForwardR = testResultForwardR, 
       testResultBackwardR = testResultBackwardR)
  
}

##
# We compute the some required values to make a decision
#

determine_causal_direction <- function(x, y) {
  x <- scale(x)
  y <- scale(y)
  
  x <- matrix(x, length(x), 1)
  y <- matrix(y, length(y), 1)
  
  result <- method_computations(x, y)
  
  # We look for the direction with the largest difference in the gaussianity of the residuals
  
  if (abs(result$testResultForward - result$testResultBackward) > 
      abs(result$testResultForwardR - result$testResultBackwardR)) {
    if (result$testResultForward > result$testResultBackward)
      list(relation = "x->y", conf = abs(result$testResultForward - result$testResultBackward))
    else
      list(relation = "y->x", conf = abs(result$testResultForward - result$testResultBackward))
  } else {
    if (result$testResultForwardR > result$testResultBackwardR)
      list(relation = "y->x", conf = abs(result$testResultForwardR - result$testResultBackwardR))
    else
      list(relation = "x->y", conf = abs(result$testResultForwardR - result$testResultBackwardR))
  }
  
}




