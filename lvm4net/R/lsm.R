#' Latent Space Model
#'
#' Latent space models (LSM) are a well known family of latent variable models for network data introduced by Hoff et al. (2002) under the basic assumption that each node has an unknown position in a D-dimensional Euclidean latent space: generally the smaller the distance between two nodes in the latent space, the greater the probability of them being connected. Unfortunately, the posterior distribution of the LSM cannot be computed analytically. For this reason we propose a variational inferential approach which proves to be less computationally intensive than the MCMC procedure proposed in Hoff et al. (2002) (implemented in the \code{latentnet} package) and can therefore easily handle large networks.
#' Salter-Townshend and Murphy (2013) applied variational methods to fit the LSM with the Euclidean distance in the \code{VBLPCM} package.
#' In this package, a distance model with squared Euclidean distance is used. We follow the notation of Gollini and Murphy (2014).
#'
#' @param Y (\code{N} x \code{N}) binary adjacency matrix
#' @param D integer dimension of the latent space
#' @param sigma (\code{D} x \code{D}) variance/covariance matrix of the prior distribution for the latent positions. Default \code{sigma = 1}
#' @param xi mean of the prior distribution of \eqn{\alpha}. Default \code{xi = 0}
#' @param psi2 variance of the prior distribution of \eqn{\alpha}. Default \code{psi2 = 2}
#' @param Niter maximum number of iterations. Default \code{Niter = 100}
#' @param Miniter minimum number of iterations. Default \code{Miniter = 10}
#' @param tol desired tolerance. Default \code{tol = 0.1^2}
#' @param randomZ logical; If \code{randomZ = TRUE} random initialization for the latent positions is used. If \code{randomZ = FALSE} and \code{D} = 2 or 3 the latent positions are initialized using the Fruchterman-Reingold method and multidimensional scaling is used for \code{D} = 1 or \code{D} > 3. Default \code{randomZ = FALSE}
#' @param nstart number of starts
#'  @return List containing:
#' \itemize{
#' \item \code{lsmEZ} (\code{N} x \code{D}) matrix containing the posterior means of the latent positions
#' \item \code{lsmVZ} (\code{D} x \code{D}) matrix containing the posterior variance of the latent positions
#' \item \code{xiT} mean of the posterior distribution of \eqn{\alpha}
#' \item \code{psi2T} variance of the posterior distribution of \eqn{\alpha}
#' \item \code{Ell} expected log-likelihood
#' }

#' @seealso \code{\link{plot.lsm}}
#' @references Gollini, I., and Murphy, T. B. (2014), "Joint Modelling of Multiple Network Views", Journal of Computational and Graphical Statistics \url{http://arxiv.org/abs/1301.3759}.
#' @references Hoff, P., Raftery, A., and Handcock, M. (2002), "Latent Space Approaches to Social Network Analysis", Journal of the American Statistical Association, 97, 1090--1098.
#' @export
#' @examples
#' ### Simulate Undirected Network
#' N <- 20
#' Y <- network(N, directed = FALSE)[,]
#'
#' modLSM <- lsm(Y, D = 2) 
#' plot(modLSM, Y)

lsm <- function(Y, D, sigma = 1, xi = 0, psi2 = 2, Niter = 100, Miniter = 10, tol = 0.1^2, randomZ = FALSE, nstart = 1)
{	
	
	Y <- as.matrix(Y)

	stopifnot(is.adjacency(Y))
	stopifnot(length(D) == 1, D > 0, D == floor(D))
	stopifnot(sigma > 0)
	stopifnot(length(xi) == 1, length(psi2) == 1, psi2 > 0)
	stopifnot(length(Miniter) == 1, Miniter > 0, Miniter == floor(Miniter), Niter > Miniter, Niter == floor(Niter))
	stopifnot(length(tol) == 1, tol > 0)
	stopifnot(is.logical(randomZ))
	stopifnot(length(nstart) == 1, nstart > 0, nstart == floor(nstart))

	best.exll <- 1
	N<-nrow(Y)
	
	miss <- anyNA(Y)
	
	if(miss){
		Yo <- Y
		miss.edges <- which(is.na(Y), TRUE)
		nmiss <- nrow(miss.edges)
		
		if(all(t(Y) == Y, na.rm = TRUE)){ # indirected
			
			Y[miss.edges[1:(nmiss / 2),]] <- rbinom(nmiss / 2, 1, .5)
			
			ind <- upper.tri(Y)
			Y[ind] <- t(Y)[ind]
						
		} else {
			
			Y[miss.edges] <- rbinom(nmiss, 1, .5)
		}
	
	}
	
	for(ns in 1:nstart) {
	
	lsm <- NULL
	
	if(randomZ){
		
		lsm$lsmEZ <- matrix(rnorm(N * D), ncol = D)
		lsm$lsmVZ <- diag(D)
		
		} else {
			
			if(D %in% 2:3){ # Fruchterman-Reingold
				
				lsm$lsmEZ <- layout.fruchterman.reingold(graph.adjacency(Y), dim = D)
				lsm$lsmEZ <- lsm$lsmEZ / apply(lsm$lsmEZ, 2, sd)
				lsm$lsmVZ <- diag(D)

			} else { # Multidimensional Scaling
				lsm$lsmEZ<-cmdscale(as.dist(1-Y), D)
				lsm$lsmVZ<-diag(D)
		}
	
	}
	
	lsm$xiT <- glm(c(Y)~c(as.matrix(dist(lsm$lsmEZ)^2)))$coeff[1]
	names(lsm$xiT) <- NULL
	lsm$psi2T <- psi2
	
	############

	iter<-0
	dif<-1
	
	l<-seq(0,0,length.out=3)

	maxtry <- 5 # try in case of cholesky decomposition error
	ntry <- 1
	
	while (iter<Niter & dif>tol)
	{
		iter<-iter+1
	
		lsm <- mainLSM(lsm$psi2T, lsm$xiT, lsm$lsmEZ, lsm$lsmVZ, Y, xi, psi2, sigma^2)
		
		## restart in case of cholesky decomposition error
		# restart up to maxtry times
		
		SI4SigmaT <- solve(diag(D) + 4 * lsm$lsmVZ)

		if(det(SI4SigmaT) < 0 & ntry < maxtry) {
			
			lsm <- restartlsm(N, D, randomZ, Y, psi2)
			
			lsm <- mainLSM(lsm$psi2T, lsm$xiT, lsm$lsmEZ, lsm$lsmVZ, Y, xi, psi2, sigma^2)
			
			iter <- 1
			diff <- 1
			l<-seq(0,0,length.out=3)
			
			ntry <- ntry + 1
		}
		
		##

		exll <- Ell(lsm$psi2T, lsm$xiT, lsm$lsmVZ, lsm$lsmEZ, Y)	
		l <- c(l[-1], exll)
		if(iter>Miniter) dif <- l[3] - l[2]
		if ( dif< -1) dif <- abs(dif) + 1
		
		if(miss){ # check diretto indirerttp
		
		if(all(t(Y) == Y)){ # indirected
			
			if(nmiss == 2){
				pmiss <- 1 / (1 + exp(- lsm$xiT + sum((lsm$lsmEZ[miss.edges[1:(nmiss / 2), 1], ] - lsm$lsmEZ[miss.edges[1:(nmiss / 2), 2], ])^2))) 	
			} else {
				pmiss <- 1 / (1 + exp(- lsm$xiT + rowSums((lsm$lsmEZ[miss.edges[1:(nmiss / 2), 1], ] - lsm$lsmEZ[miss.edges[1:(nmiss / 2), 2], ])^2))) 	
				}
			
			Y[miss.edges[1:(nmiss / 2),]] <- pmiss
			
			ind <- upper.tri(Y)
			Y[ind] <- t(Y)[ind]
						
		} else {
			
			if(nmiss == 1){
				pmiss <- 1 / (1 + exp(- lsm$xiT + sum((lsm$lsmEZ[miss.edges[,1], ] - lsm$lsmEZ[miss.edges[,2], ])^2))) 
			} else {
				pmiss <- 1 / (1 + exp(- lsm$xiT + rowSums((lsm$lsmEZ[miss.edges[,1], ] - lsm$lsmEZ[miss.edges[,2], ])^2))) 
				}
			
			Y[miss.edges] <- pmiss
		}
				
		}
	
	}
	
	if(is.null(rownames(lsm$lsmEZ))) rownames(lsm$lsmEZ) <- paste("Obs_", seq(1:nrow(lsm$lsmEZ)), sep="")
	
	if(ns == 1 || exll > best.exll){
		best.lsm <- lsm
		best.exll <- exll
	}
	
	}
		robj <- list(lsmEZ = best.lsm$lsmEZ, lsmVZ = best.lsm$lsmVZ, xiT = best.lsm$xiT, psi2T = best.lsm$psi2T, Ell = best.exll)
	
	class(robj) <- c("lsm")
	
	robj
}