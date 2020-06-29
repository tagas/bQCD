library(energy)

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


gaussian_kernel_rkhse <- function(x, y, sigma, lambda) {

	L <- nrow(x)
	
	ux <- matrix(runif(L), L, 1)
	uy <- matrix(runif(L), L, 1)

	Kx <- gaussian_kernel(x, sigma)
	Kux <- gaussian_kernel(ux, sigma)
	Kuy <- gaussian_kernel(uy, sigma)
	Ky <- gaussian_kernel(y, sigma)
	Kxux <- gaussian_kernel_XY(x, ux, sigma)
	Kyuy <- gaussian_kernel_XY(y, uy, sigma)

	Vx <- chol2inv(chol(Kx + diag(nrow(Kx)) * lambda * L))
	Vy <- chol2inv(chol(Ky + diag(nrow(Ky)) * lambda * L))

	first <- 1 / L^2 * sum(diag(Kx %*% Ky)) - 2 / L^2 * sum(diag(Ky %*% Vx %*% Kxux %*% t(Kxux))) + 
		1 / L^2 * sum(diag(Ky %*% Vx %*% Kxux %*% Kux %*% t(Kxux) %*% Vx))
	second <- 1 / L^2 * sum(diag(Ky %*% Vx %*% Kxux %*% Kux %*% t(Kxux) %*% Vx)) - 2 / L^2 * 
		sum(diag(t(Kyuy) %*% Vx %*% Kxux %*% Kux)) + 1 / L^2 * sum(diag(Kux %*% Kuy))

	third <- 1 / L^2 * sum(diag(Ky %*% Kx)) - 2 / L^2 * sum(diag(Kx %*% Vy %*% Kyuy %*% t(Kyuy))) + 
		1 / L^2 * sum(diag(Kx %*% Vy %*% Kyuy %*% Kuy %*% t(Kyuy) %*% Vy))
	fourth <- 1 / L^2 * sum(diag(Kx %*% Vy %*% Kyuy %*% Kuy %*% t(Kyuy) %*% Vy)) - 2 / L^2 * 
		sum(diag(t(Kxux) %*% Vy %*% Kyuy %*% Kuy)) + 1 / L^2 * sum(diag(Kuy %*% Kux))

	Cxy <- first + second
	Cyx <- third + fourth

	list(Cxy = Cxy, Cyx = Cyx)

}

linear_kernel_rkhse <- function(x, y, lambda) {

	L <- nrow(x)
	
	ux <- matrix(runif(L), L, 1)
	uy <- matrix(runif(L), L, 1)

	Kx <- linear_kernel(x)
	Kux <- linear_kernel(ux)
	Kuy <- linear_kernel(uy)
	Ky <- linear_kernel(y)
	Kxux <- linear_kernel_XY(x, ux)
	Kyuy <- linear_kernel_XY(y, uy)

	Vx <- chol2inv(chol(Kx + diag(nrow(Kx)) * lambda * L))
	Vy <- chol2inv(chol(Ky + diag(nrow(Ky)) * lambda * L))

	first <- 1 / L^2 * sum(diag(Kx %*% Ky)) - 2 / L^2 * sum(diag(Ky %*% Vx %*% Kxux %*% t(Kxux))) + 
		1 / L^2 * sum(diag(Ky %*% Vx %*% Kxux %*% Kux %*% t(Kxux) %*% Vx))
	second <- 1 / L^2 * sum(diag(Ky %*% Vx %*% Kxux %*% Kux %*% t(Kxux) %*% Vx)) - 2 / L^2 * 
		sum(diag(t(Kyuy) %*% Vx %*% Kxux %*% Kux)) + 1 / L^2 * sum(diag(Kux %*% Kuy))

	third <- 1 / L^2 * sum(diag(Ky %*% Kx)) - 2 / L^2 * sum(diag(Kx %*% Vy %*% Kyuy %*% t(Kyuy))) + 
		1 / L^2 * sum(diag(Kx %*% Vy %*% Kyuy %*% Kuy %*% t(Kyuy) %*% Vy))
	fourth <- 1 / L^2 * sum(diag(Kx %*% Vy %*% Kyuy %*% Kuy %*% t(Kyuy) %*% Vy)) - 2 / L^2 * 
		sum(diag(t(Kxux) %*% Vy %*% Kyuy %*% Kuy)) + 1 / L^2 * sum(diag(Kuy %*% Kux))

	Cxy <- first + second
	Cyx <- third + fourth

	list(Cxy = Cxy, Cyx = Cyx)
}

find_sigma_median <- function(x, y) {

	distances_x <- as.matrix(dist(x))
	distances_x[ lower.tri(distances_x) ] <- 0
	distances_y <- as.matrix(dist(y))
	distances_y[ lower.tri(distances_y) ] <- 0

	sigma <- 1 / (2 * median(c(distances_x[ distances_x > 0 ], distances_y[ distances_y > 0 ]))) 

	sigma * 5
}

determine_causal_direction_rkhse <- function(x, y) {

	x <- (x - min(x)) / (max(x) - min(x))
	y <- (y - min(y)) / (max(y) - min(y))
	
	x <- matrix(x, length(x), 1)
	y <- matrix(y, length(y), 1)

	sigma <- find_sigma_median(x, y)
	lambda <- 1e-3

	result <- gaussian_kernel_rkhse(x, y, sigma, lambda)

	if (result$Cxy < result$Cyx)
		list(relation = "x->y", conf = abs(result$Cxy - result$Cyx))
	else
		list(relation = "y->x", conf = abs(result$Cxy - result$Cyx))
}




