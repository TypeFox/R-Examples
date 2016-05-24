#' This function returns the synthetic spatiotemporal data set resembling functional MR Images (fMRI) data.
#'
#' The returned data is simulated on a 20 x 20 grid. 
#'
#' @name sim.fmri2COVAR
#' @aliases sim.fmri2COVAR
#' @title Simulate FMRI Data
#' @usage sim.fmri2COVAR(hrf, beta.Var1, beta.Var2)
#' @param hrf haemodynamic response function, needs to be a vector of length \code{T}.
#' @param beta.Var1 scalar, defines the height of the activated area, in form of a cylinder of the first grid. 
#' @param beta.Var2 scalar, defines the height of the activated area, in form of a cylinder of the second grid. 
#' @author Max Hughes
#' @note This function is solely for two covariates.
#' @examples
#' # non-transformed hr-function
#' T <- 180
#' seq.length <- T*3
#' index <- seq(3, T*3, by = 3)
#' vis <- rep(c(-0.5, 0.5), each=30, times=ceiling(T/30*1.5))
#' vis <- as.matrix(vis[index])
#' aud <- rep(c(-0.5, 0.5), each=45, times=ceiling(T/30*1.5))
#' aud <- as.matrix(aud[index])
#' hrf <- cbind(vis,aud)
#' # define height of activation area
#' beta.Var1 <- beta.Var2 <- 3
#' # use function to obtain fmri data
#' data <- sim.fmri2COVAR(hrf, beta.Var1, beta.Var2)$fmri            


sim.fmri2COVAR <- function(hrf, beta.Var1, beta.Var2){

  require(Matrix)

  if(any(is.na(hrf)))
    stop("\nNAs in hr function.\n")

  ## For model with two covariates
  p <- dim(hrf)[2]
  if(p!=2)
    stop("Haemodynamic response function needs to be a matrix with column dimension of 2.")
    
  I <- 400
  ## covariate(s) (for now only one: visual)
  Z <- as.matrix(hrf)        # note: Z is only for one pixel
  T <- dim(Z)[1]
  Z.Var1 <- as.matrix(Z[,1])
  Z.Var2 <- as.matrix(Z[,2])
  IZ.Var1 <- kronecker(as(diag(1, nrow=I, ncol=I), "sparseMatrix"), Z.Var1)
  IZ.Var2 <- kronecker(as(diag(1, nrow=I, ncol=I), "sparseMatrix"), Z.Var2)

  ## generate sigma^2
  sigma.sq <- numeric(I)
  for(i in 1:I)
    sigma.sq[i] <- 25 + 2*rnorm(1)

  ## generate epsilon
  eps <- matrix(nrow=T, ncol=I)
    for(i in 1:I)
      eps[,i] <- rnorm(T, mean=0, sd=sqrt(sigma.sq[i]))    # Variance or sd???

  eps <- as.vector(eps)

  ## generate beta grid for Variable 1
  beta.sim.Var1 <- matrix(0, nrow=I/20, ncol=I/20)
  beta.sim.Var1[8:14,8:14] <- beta.Var1

  beta.sim.Var1[8,8] <- 0
  beta.sim.Var1[9,8] <- 0
  beta.sim.Var1[8,9] <- 0

  beta.sim.Var1[13,8] <- 0
  beta.sim.Var1[14,9] <- 0
  beta.sim.Var1[14,8] <- 0

  beta.sim.Var1[9,14] <- 0
  beta.sim.Var1[8,13] <- 0
  beta.sim.Var1[8,14] <- 0

  beta.sim.Var1[14,13] <- 0
  beta.sim.Var1[13,14] <- 0
  beta.sim.Var1[14,14] <- 0

  beta.sim.Var1 <- as.vector(beta.sim.Var1)
  
  ## generate beta grid for Variable 2
  beta.sim.Var2 <- matrix(0, nrow=I/20, ncol=I/20)
  beta.sim.Var2[8:14,8:14] <- beta.Var2

  beta.sim.Var2[8,8] <- 0
  beta.sim.Var2[9,8] <- 0
  beta.sim.Var2[8,9] <- 0

  beta.sim.Var2[13,8] <- 0
  beta.sim.Var2[14,9] <- 0
  beta.sim.Var2[14,8] <- 0

  beta.sim.Var2[9,14] <- 0
  beta.sim.Var2[8,13] <- 0
  beta.sim.Var2[8,14] <- 0

  beta.sim.Var2[14,13] <- 0
  beta.sim.Var2[13,14] <- 0
  beta.sim.Var2[14,14] <- 0

  beta.sim.Var2 <- as.vector(beta.sim.Var2)

  ## generate y
  y <- numeric(T*I)
  y <- IZ.Var1%*%beta.sim.Var1+IZ.Var2%*%beta.sim.Var2+eps

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

  return(list("fmri"=y, "hrf"=Z, "coeff1"=beta.sim.Var1, "coeff2"=beta.sim.Var2,
              "resid"=eps, "sigma"=sigma.sq))
}

