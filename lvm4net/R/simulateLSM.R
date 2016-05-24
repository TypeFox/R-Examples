#' Simulate from LSM model
#'
#' Function to simulate networks from the LSM model
#'
#' @param object object of class \code{'lsm'}
#' @param Y (\code{N} x \code{N}) binary adjacency matrix
#' @param nsim number of simulations. Default \code{nsim = 100}
#' @param seed for simulations 
#' @param directed if the network is directed or not.  Default \code{directed = NULL} 
#' @export
#' @examples
#' 
#' n <- 20
#' Y <- network(n, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#'
#' Ysim <- simulateLSM(modLSM, Y = Y, nsim = 8)
#' # store EZ, to keep the nodes in the same positions 
#' # and compare the networks
#' EZ <- modLSM$lsmEZ
#' par(mfrow = c(3,3))
#' plotY(Y, EZ = EZ, main = "Original Data")
#' for(i in 1:8) plotY(Ysim[[i]], EZ = EZ, main = paste("Simulation" , i))
#' par(mfrow = c(1,1))

simulateLSM <- function(object, Y = NULL, nsim = 100, seed, directed = NULL){
	
	stopifnot(inherits(object, 'lsm'))
	
	nsim <- as.integer(nsim)
	stopifnot(nsim > 0)
	
	xiT <- object$xiT
	psiT <- sqrt(object$psi2T)
	
	EZ <- object$lsmEZ
	VZ <- object$lsmVZ
	
	N <- nrow(EZ)
	
	if(is.null(directed)){ 
			directed <- !all(Y == t(Y))
		} else {
			stopifnot(is.logical(directed) & length(directed) == 1)
		}
	
	Ysim <- list("Ysim")
	
	if(directed){
		
	for (i in 1:nsim) {
		
		xiTs <- rnorm(1, mean = xiT, sd = psiT)
		EZs <- apply(EZ, 1, function(ez) mvrnorm(mu = ez, Sigma = VZ))
		num <- exp(xiTs - dist(t(EZs))^2)
		p1 <- as.matrix(num / (1 + num))
		
		Ysim[[i]] <- matrix(rbinom(N * N, 1, p1), N, N)
		diag(Ysim[[i]]) <- 0
	}
	
	} else {
	 
	 for (i in 1:nsim){
	 	
	 	xiTs <- rnorm(1, mean = xiT, sd = psiT)
		EZs <- apply(EZ, 1, function(ez) mvrnorm(mu = ez, Sigma = VZ))
		num <- exp(xiTs - dist(t(EZs))^2)
		p1 <- c(num / (1 + num))
	 	
		y <- matrix(0, N, N)
		y[lower.tri(y)] <- rbinom(N * (N - 1) / 2, 1, prob = p1)
		y[upper.tri(y)] <- t(y)[upper.tri(y)]
		
		Ysim[[i]] <- y 
		}
	
	}
	
	Ysim
}