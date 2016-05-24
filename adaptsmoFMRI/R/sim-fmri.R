#' This function returns the synthetic spatiotemporal data set resembling functional MR Images (fMRI) data.
#'
#' The returned data is simulated on a 20 x 20 grid. 
#'
#' @name sim.fmri 
#' @aliases sim.fmri
#' @title Simulate FMRI Data
#' @usage sim.fmri(hrf, beta)
#' @param hrf haemodynamic response function, needs to be a vector of length \code{T}.
#' @param beta scalar, defines the height of the activated area, in form of a cylinder. 
#' @author Max Hughes
#' @note This function is solely for one covariate.
#' @examples
#' # non-transformed hr-function
#' T <- 210
#' seq.length <- T*3
#' index <- seq(3, T*3, by = 3)
#' hrf <- rep(c(-0.5, 0.5), each=30, times=ceiling(T/30*1.5))
#' hrf <- as.matrix(hrf[index])
#' # define height of activation area
#' beta <- 3
#' # use function to obtain fmri data
#' data <- sim.fmri(hrf, beta)$fmri                   

sim.fmri <- function(hrf, beta){

  require(Matrix)

  if(any(is.na(hrf)))
    stop("\nNAs in hr function.\n")
        
  I <- 400
  ## covariate(s) (for now only one: visual)
  Z <- as.matrix(hrf)        # note: Z is only for one pixel
  T <- dim(Z)[1]
  IZ <- kronecker(as(diag(1, nrow=I, ncol=I), "sparseMatrix"), Z)

  ## generate sigma^2
  sigma.sq <- numeric(I)
  for(i in 1:I)
    sigma.sq[i] <- 25 + 2*rnorm(1)

  ## generate epsilon
  eps <- matrix(nrow=T, ncol=I)
    for(i in 1:I)
      eps[,i] <- rnorm(T, mean=0, sd=sqrt(sigma.sq[i]))    # Variance or sd???

  eps <- as.vector(eps)

  ## generate beta grid
  beta.sim <- matrix(0, nrow=I/20, ncol=I/20)
  beta.sim[8:14,8:14] <- beta
  
  beta.sim[8,8] <- 0
  beta.sim[9,8] <- 0
  beta.sim[8,9] <- 0

  beta.sim[13,8] <- 0
  beta.sim[14,9] <- 0
  beta.sim[14,8] <- 0

  beta.sim[9,14] <- 0
  beta.sim[8,13] <- 0
  beta.sim[8,14] <- 0

  beta.sim[14,13] <- 0
  beta.sim[13,14] <- 0
  beta.sim[14,14] <- 0

  beta.sim <- as.vector(beta.sim)
 
  ## generate y
  y <- numeric(T*I)
  y <- IZ%*%beta.sim+eps

  ## convert back to to an array of dim (nrow x ncol x T), so that it coincides with
  ## the dim of real data
  y <- t(matrix(nrow=T, ncol=I, data=y))
  y <- array(y, dim=c(20,20,T))
  # interchange the first two subscripts on a 3-way array y
  y <- aperm(y, c(2,1,3))

  # same procededure for epsilon
  eps <- t(matrix(nrow=T, ncol=I, data=eps))
  eps <- array(eps, dim=c(20,20,T))
  eps <- aperm(eps, c(2,1,3))
  
  return(list("fmri"=y, "hrf"=Z, "coeff"=beta.sim, "resid"=eps,
              "sigma"=sigma.sq))
}




