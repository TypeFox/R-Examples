
# (C) 2008-2012 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "Where a calculator on the ENIAC is equipped with 18,000 vacuum tubes
# and weighs 30 tons, computers in the future may have only 1,000
# vaccuum tubes and perhaps weigh 1.5 tons."
# - unknown, Popular Mechanics, March 1949


centerData <- function (X, rm.na = TRUE, meanvalue = NULL) {

  # Shift data matrix (columns) to zero, or given 'meanvalue'
  
  if (!rm.na) {
    xcenter <- colMeans(X)
    X2 <- X - rep(xcenter, rep.int(nrow(X), ncol(X)))
  } else {	
    X2 <- array(NA, dim = c(nrow(X), ncol(X)), dimnames = dimnames(X))
    for (i in 1:ncol(X)) {
      x <- X[,i]
      nainds <- is.na(x)
      xmean <- mean(x[!nainds])
      X2[!nainds,i] <- x[!nainds] - xmean 	
    }
    dimnames(X2) <- dimnames(X)
  }

  if (!is.null(meanvalue)) {
    # Shift the data so that mean gets a specified value
    X2 <- X2 + meanvalue
  }

  X2
}


initialize2 <- function (X, Y, zDim = NULL, marginalCovariances) {

  zDim <- ifelse(is.null(zDim), min(nrow(X), nrow(Y)), zDim)

  Nsamples <- ncol(X)
  Dim <- list(X = nrow(X), Y = nrow(Y), Z = zDim)
  nullmat  <- matrix(0, nrow = Dim$X, ncol = Dim$Y)

  if (marginalCovariances == "isotropic") {
    # Scalar values
    phi <- list()
    phi$X <- diag(var(as.vector(X)), Dim$X)
    phi$Y <- diag(var(as.vector(Y)), Dim$Y)
    phi$total <- diag(c(diag(phi$X), diag(phi$Y)))
    #nullmat <- 0
  } else if (marginalCovariances == "identical isotropic") {
    # Scalar values phix = phiy
    phi <- list()
    phi.est <- var(c(as.vector(X), as.vector(Y)))
    phi$X <- diag(phi.est, Dim$X)
    phi$Y <- diag(phi.est, Dim$Y)
    phi$total <- diag(c(diag(phi$X), diag(phi$Y)))    
    #nullmat <- 0
  } else {
    # diagonal matrices
    # initialize with scalar diagonal noise on the marginals (shared by all features)
    phi <- list(X = diag(var(as.vector(X)), Dim$X), 
                Y = diag(var(as.vector(Y)), Dim$Y))  
    phi$total <- rbind(cbind(phi$X,nullmat), cbind(t(nullmat), phi$Y))

  }

  # FIXME: if phi$Y is scalar (as in segmented/mir case) we can speed up here. Do later.

  phi.inv  <- list()
  phi.inv$X <- solve(phi$X)
  phi.inv$Y <- solve(phi$Y)
  #phi.inv$total <- rbind(cbind(phi.inv$X, nullmat), cbind(t(nullmat), phi.inv$Y))

  Dcov <- list()
  Dcov$X <- cov(t(X), use = "pairwise.complete.obs")
  Dcov$Y <- cov(t(Y), use = "pairwise.complete.obs")
  Dcov$xy <- cov(t(X), t(Y), use = "pairwise.complete.obs")
  Dcov$yx <- t(Dcov$xy)
  Dcov$total <- rbind(cbind(Dcov$X, Dcov$xy), cbind(Dcov$yx, Dcov$Y))
  if (nrow(X) == nrow(Y)) {
    Dcov$sum <- Dcov$X + Dcov$Y + Dcov$xy + Dcov$yx
    #Dcov$sum <- cov(t(X + Y), use = "pairwise.complete.obs")
  }

  # It is possible that covariances calculated with pairwise complete
  # observations are not positive semi-definite.
  # Check this. If not pos.sem.def, then replace with the closest
  # pos. semidefinite matrix
  if (any(eigen(Dcov$total)$values < 0)) {
    #message("Covariance approximation used to avoid numerical instability.")
    Dcov$X  <- as.matrix(nearPD(Dcov$X)$mat)
    Dcov$Y  <- as.matrix(nearPD(Dcov$Y)$mat)
    Dcov$total <- rbind(cbind(Dcov$X, Dcov$xy), cbind(Dcov$yx, Dcov$Y))
    if (nrow(X) == nrow(Y)) {
      Dcov$sum <- Dcov$X + Dcov$Y + Dcov$xy + Dcov$yx
    }
  }

  # Initialize W's
  W <- list()
  W$X <- as.matrix(eigen(Dcov$X)$vectors[, 1:Dim$Z])
  W$Y <- as.matrix(eigen(Dcov$Y)$vectors[, 1:Dim$Z])	
  W$total <- rbind(W$X, W$Y) 
  
  list(phi = phi, W = W, Dcov = Dcov, Dim = Dim, nullmat = nullmat, Nsamples = Nsamples)

}

