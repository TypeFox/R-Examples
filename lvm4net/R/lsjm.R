#' Latent Space Joint Model
#'
#' Function to joint modelling of multiple network views using the Latent Space Jont Model (LSJM) Gollini and Murphy (2014). 
#' The LSJM merges the information given by the multiple network views by assuming that the probability of a node being connected with other nodes in each view is explained by a unique latent variable.
#'
#' @param Y list containing a (\code{N} x \code{N}) binary adjacency matrix for each network view.
#' @param D integer dimension of the latent space
#' @param sigma (\code{D} x \code{D}) variance/covariance matrix of the prior distribution for the latent positions. Default \code{sigma = 1}
#' @param xi vector of means of the prior distributions of \eqn{\alpha}. Default \code{xi = 0}
#' @param psi2 vector of variances of the prior distributions of \eqn{\alpha}. Default \code{psi2 = 2}
#' @param Niter maximum number of iterations. Default \code{Niter = 500}
#' @param tol desired tolerance. Default \code{tol = 0.1^2}
#' @param preit Preliminary number of iterations default \code{preit = 20}
#' @param randomZ logical; If \code{randomZ = TRUE} random initialization for the latent positions is used. If \code{randomZ = FALSE} and \code{D} = 2 or 3 the latent positions are initialized using the Fruchterman-Reingold method and multidimensional scaling is used for \code{D} = 1 or \code{D} > 3. Default \code{randomZ = FALSE}

#'  @return List containing:
#' \itemize{
#' \item \code{EZ} (\code{N} x \code{D}) matrix containing the posterior means of the latent positions
#' \item \code{VZ} (\code{D} x \code{D}) matrix containing the posterior variance of the latent positions
#' \item \code{lsmEZ} list contatining a (\code{N} x \code{D}) matrix for each network view containing the posterior means of the latent positions under each model in the latent space.
#' \item \code{lsmVZ} list contatining a (\code{D} x \code{D}) matrix for each network view containing the posterior variance of the latent positions under each model in the latent space.
#' \item \code{xiT} vector of means of the posterior distributions of \eqn{\alpha}
#' \item \code{psi2T} vector of variances of the posterior distributions of \eqn{\alpha}
#' \item \code{Ell} expected log-likelihood
#' }
#' @references Gollini, I., and Murphy, T. B. (2014), "Joint Modelling of Multiple Network Views", Journal of Computational and Graphical Statistics \url{http://arxiv.org/abs/1301.3759}.
#' @export
#' @examples
#'## Simulate Undirected Network
#'   N <- 20
#'   Ndata <- 2
#'    Y <- list()
#'    Y[[1]] <- network(N, directed = FALSE)[,]
#'    ### create a new view that is similar to the original
#'    
#'   for(nd in 2:Ndata){
#'     Y[[nd]] <- Y[[nd - 1]] - sample(c(-1, 0, 1), N * N, replace = TRUE, 
#'    prob = c(.05, .85, .1))
#'     Y[[nd]] <- 1 * (Y[[nd]]  > 0 )
#'   diag(Y[[nd]]) <- 0
#'    }
#'
#' par(mfrow = c(1, 2))
#' z <- plotY(Y[[1]], verbose = TRUE, main = 'Network 1')
#' plotY(Y[[2]], EZ = z, main = 'Network 2')
#' par(mfrow = c(1, 1))
#'
#' modLSJM <- lsjm(Y, D = 2) 
#' plot(modLSJM, Y, drawCB = TRUE)
#' plot(modLSJM, Y, drawCB = TRUE, plotZtilde = TRUE)

lsjm<-function(Y, D, sigma = 1, xi = rep(0, length(Y)), psi2 = rep(2, length(Y)), 
               Niter = 500, tol = 0.1^2, preit = 20, randomZ = FALSE)
{
	stopifnot(is.list(Y), sapply(Y, is.adjacency))
	stopifnot(length(D) == 1, D > 0, D == floor(D))
	stopifnot(sigma > 0)
	stopifnot(length(xi) == length(Y), length(psi2) == length(Y), psi2 > 0)
	stopifnot(preit > 0, preit == floor(preit), Niter > preit, Niter == floor(Niter))
	stopifnot(tol > 0)
	stopifnot(is.logical(randomZ))
	
	stepA <- 0
	N <- nrow(Y[[1]])
	Ndata <- length(Y)
	
	xiT <- xi
	psi2T<-psi2
	
	lsmVZ <- list()
	
	for(i in 1:Ndata) lsmVZ[[i]] <- diag(D)

	if(randomZ){
		
		lsmEZ <- list()
		for(i in 1:Ndata) lsmEZ[[i]]  <- matrix(rnorm(N * D), ncol = D)
		
		} else {
			
			if(D %in% 2:3){ # Fruchterman-Reingold
				lsmEZ <- lapply(Y, frEZ, d = D)
			} else { # Multidimensional Scaling
				lsmEZ <- lapply(Y, function(y) cmdscale(as.dist(1 - y), k = D))
		}
	
	}
	
	lsmEZ <- lapply(lsmEZ, function(z) z / apply(z, 2, sd))
	
	Aezvz <- 0

	for(i in 1:Ndata)
	{		
		Aezvz <- lsmEZ[[i]] %*% solve(lsmVZ[[i]]) + Aezvz
		
		xiT[i] <- glm(c(Y[[i]])~c(as.matrix(dist(lsmEZ[[i]])^2)))$coeff[1]
		names(xiT[i]) <- NULL
	}
	
if(D>1){
	
		VZ <- solve(matrix(rowSums(sapply(lsmVZ, solve)), D, D) - (Ndata - 1) / sigma^2 * diag(D))
	
	} else {
		VZ <- as.matrix(1 / sum(sapply(lsmVZ, solve)) - (Ndata - 1) / sigma^2)
	}
	
	EZ <- Aezvz %*% VZ	
	
	############
	############

	iter <- 0
	dif <- 1
	l <- seq(0, 0, length.out=3)
	ellm <- rep(0, Ndata)

	while (iter<Niter & dif>tol)
	{
		iter<- iter+1	
	
	for(i in 1:Ndata)
	{	
		lsm <- mainLSM(psi2T[i], xiT[i], EZ, VZ, Y[[i]], xi[i], psi2[i], sigma^2)
		xiT[i] <- lsm$xiT
		psi2T[i] <- lsm$psi2T
		lsmVZ[[i]] <- as.matrix(lsm$lsmVZ)
		lsmEZ[[i]] <- lsm$lsmEZ
	}

	######## Joint Model ##############	

	if(D > 1 & iter > preit)
	{
	for(i in 2:Ndata)
 	{	
 		rotLSM <- rotXtoY(lsmEZ[[i]], lsmEZ[[1]] %*% lsmVZ[[i]])
 		lsmEZ[[i]] <- rotLSM$X
 	}
	}
	
	Aezvz <- 0
	
	for(i in 1:Ndata) {	
 		Aezvz <- lsmEZ[[i]] %*% solve(lsmVZ[[i]]) + Aezvz
 	}
 	
 	#####

	if(D>1){
		VZ <- solve(matrix(rowSums(sapply(lsmVZ, solve)), D, D) - (Ndata - 1) / sigma^2 * diag(D))
	} else {
		VZ <- as.matrix(1 / sum(sapply(lsmVZ, solve)) - (Ndata - 1) / sigma^2)
	}
	
	EZ<- Aezvz %*% VZ

	########
	
	for(i in 1:Ndata)
	{	
		ellm[i] <- Ell(psi2T[i], xiT[i], VZ, EZ, Y[[i]])
	}

	ell <-sum(ellm)	
	l <- c(l[-1], ell)
	
	if(iter > preit) dif<- l[3] - l[2]
		
	if(dif < -tol) dif <- abs(dif) + 1
	
	}
	
	robj <- list(EZ = EZ, VZ = VZ, lsmEZ = lsmEZ, lsmVZ = lsmVZ, xiT = xiT, psi2T = psi2T, Ell = ell)
	
	class(robj) <- c("lsjm")
	
	robj
	
} 

