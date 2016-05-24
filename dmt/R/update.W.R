# (C) 2008-2012 Leo Lahti and Olli-Pekka Huovilainen          
# All rights reserved. 
# FreeBSD License (keep this notice)     


# "The mathematics is not there till we put it there."
# - Sir Arthur Eddington, The Philosophy of Physical Science


update.W.singledata <- function (Wt, X, phi, priors = NULL) {

  # Optimize W for PFA and PPCA models
  
  # Use cost function to add priors on W
  if (is.null(priors)) {
    Wvec <- optim(as.vector(Wt), pfa.neg.log.likelihood, X = X, phi = phi, method = "SANN")$par # NelderMead and BFGS do not work as well 
  } else if (!is.null(priors$W)) {
    Wvec <- abs(optim(as.vector(Wt), pfa.cost.regularized, X = X, phi = phi, priors = priors, method = "SANN")$par) # NelderMead and BFGS do not work as well 
    
  }
  
  W <- matrix(Wvec, ncol = nrow(X))

}


pfa.cost.regularized <- function (Wvec, phi, X, priors) {

  Wvec <- abs(Wvec)
  
  cost.data <- pfa.neg.log.likelihood(Wvec, phi, X) 

  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  #wcost <- 0 # no effect # FIXME: would be faster without this 'if' check  
  #if ( !is.null(priors$W) && priors$W > 0 ) {
  #wcost <- -sum(dexp(Wvec, rate = priors$W, log = TRUE))
  #wcost <- wprior.c(Wvec, priors$W)
  wcost <- wprior(Wvec, priors$W)
  #}
  
  cost.data + wcost
    
}

wprior <- function (vec, rate) {
  -sum(dexp(vec, rate = rate, log = TRUE))
}
#wprior.c <- cmpfun(wprior)


pfa.neg.log.likelihood <- function (Wvec, phi, X) {

  # Cost function for W in the PFA model. 
  # Also applicable for PPCA.
  
  W <- matrix(Wvec, ncol = nrow(X))

  # X: features x samples; assuming that this is centered at origo

  # X is Y in Rubin-Thayer 1982: this log-likelihood is from Eq. 1 in there
  # R <- diag(1, zDimension) # R = I ie. exploratory factor analysis, see Rubin-Thayer Case 1.
  # k <- tau2 + t(beta) %*% R %*% beta
  # k <- tau2 + t(beta) %*% beta # assuming R = I
  # beta <- t(W$total)
  # tau2 <- phi$total
  # k <- tau2 + t(beta) %*% beta # assuming R = I
  # wtw <- W%*%t(W)
  
  # assuming R = I and adding small constant to avoid numerical overflows
  k <- t(W)%*%W + phi + diag(1e-18, nrow(phi)) 
  pfacost(ncol(X), k, X)
  #pfacost.c(ncol(X), k, X)

}

pfacost <- function (n, k, X) {
  -as.numeric(-(n/2)*(determinant(k, logarithm = TRUE)$modulus + sum(diag(cov(t(X)) %*% solve( k )))))
}



cost.W <- function (vec, phi, priors, Dim, Dcov) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # Wx ~ Wy constrained
  # no W prior

  # Retrieve the actual W and T from the parameter vector
  wt <- get.W(vec, Dim)
  W <- wt$W
  T <- wt$T

  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # We report -logP here  
  # Data prob. Taken from probCCA paper, section 4, l1

  wtw.xy <- W$X%*%t(W$Y)
  Sigma <- rbind(cbind(W$X%*%t(W$X) + phi$X*diag(Dim$X), wtw.xy),
                 cbind(t(wtw.xy),W$Y%*%t(W$Y) + phi$Y*diag(Dim$Y)))

  # -logP for the data
  cost.data <- log(det(Sigma)) + sum(diag(solve(Sigma)%*%Dcov$total))

  # -logP for T prior
  tcost <- sum((T - priors$Nm.wxwy.mean)^2) * priors$T.tmp

  # -logP for W prior - skip since not used now
  # priors$W.tmp <- 1/(2 * Nsamples * priors$W)
  # NOTE considerable speed increase in optimize iteration if 
  # this is calculated outside this function!
  #wcost <- sum((W$X)^2) * priors$W.tmp

  cost.data + tcost #+ wcost

}

cost.W.exponential <- function (vec, phi, priors = NULL, Dim, Dcov) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # allows exponential prior for W
  # in general, Wx != Wy

  # remove sign as we assume W always positive here
  vec <- abs(vec)
  
  # Retrieve W from the parameter vector
  W <- get.W.nonneg(vec, Dim)

  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # report -logP here
  
  # Data prob. Taken from probCCA paper, section 4, l1
  wtw.xy <- W$X%*%t(W$Y)

  Sigma <- rbind(cbind(W$X%*%t(W$X) + phi$X, wtw.xy),
                 cbind(t(wtw.xy), W$Y%*%t(W$Y) + phi$Y))
  
  # -logP for the data
  det.sigma <- det(Sigma) 

  if (det.sigma > 1e-6) { # using > 0 caused overflows 
    cost.data <- log(det.sigma) + sum(diag(solve(Sigma)%*%Dcov$total))
  } else {
    # In some rare situations det.sigma appears non-positive;
    # such solutions are not feasible
    cost.data <- 1e300 # 1e309 gives Inf; Inf gives error
  }
  
  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  wcost <- 0 # no effect # FIXME: would be faster without this 'if' check
  if ( !is.null(priors$W) && priors$W > 0 ) {
    wcost <- -sum(dexp(vec, rate = priors$W, log = TRUE))
  }
  
  cost.data + wcost

}

cost7 <- function (Wvec, phi, Dcov, Dim, priors) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # SimCCA: identical Wx = Wy
  # allows W prior

  # NOTE: possible to optimize quite much by removing W matrix conversions?

  if ( !is.null(priors$W) ) { Wvec <- abs(Wvec) }    

  W <- get.W.nonneg.identical(Wvec, Dim)
  wtw <- W%*%t(W)

  Sigma <- rbind(cbind(wtw + phi$X, wtw),
                 cbind(wtw, wtw + phi$Y))


  # Marginal cost for the whole data set
  # integrated over z
  # given parameters W, phi
  # P(X,Y | W, phi) = integral N(X|Wx*z,phix)*N(Y|Wy*z,phiy)*N(z|0,I)
  # reporting -logP here

  # restrict solutions to cases where det(Sigma)>=0

  # -logP for the data
  detsigma <- det( Sigma )

  if (detsigma > 1e-100) { # crashes with > 0
    cost.data <- log(detsigma) + sum(diag(solve(Sigma)%*%Dcov$total))
  } else { cost.data <- 1e300 }
    
  # -logP for W prior
  # wcost <- sum((W$X)^2) * priors$W
  # Assuming exponential prior distribution with rate parameter priors$W
  wcost <- 0 # no effect by efault
  if ( !is.null(priors$W) ) {
    #multiply by 2 to count for both wx and wy
    wcost <- -2*sum(dexp(Wvec, rate = priors$W, log = TRUE))
  } 

  cost.data + wcost

}


W.cca.EM <- function (Dcov, M, beta) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # EM update for W, given phi (through M, beta)
  # Returns total W, i.e. [Wx; Wy]

  #beta <-M%*%t(W$total)%*%phi.inv$total
  #M <- solve(t(W)%*%phi.inv%*%W + I) # W$total meant here

  ctb <- Dcov$total%*%t(beta)
  matrix(ctb%*%solve(M + beta%*%ctb), nrow = nrow(Dcov$total))

}


W.simcca.EM <- function (W, phi, Dim, Dcov) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # CCA update for W, assuming Wx = Wy
  # see Bach-Jordan 2005, sec. 4.1 for details
  # equations modified from there to match Wx = Wy case
  # FIXME: speedup by sharnig M/beta with phi updates as in in W.cca.EM?
  # (phi.EM.simcca)
  
  what <- 2*W$X
  twp <- t(what)%*%solve(phi$X + phi$Y)
  M.w <- solve(twp%*%what + diag(Dim$Z))
  beta.w <- M.w%*%twp
  ctb <- Dcov$sum%*%t(beta.w)
  w <- ctb%*%solve(M.w + beta.w%*%ctb)/2
  
  list(X = w, Y = w, total = rbind(w, w))

}

					
#######################################
	
solve.w <- function (Xc, Yc, Cxx, Cyy, dz = NULL) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # FIXME: compare with the other W updates, e.g. W.cca.EM

  # assumes Xc, Yc : samples x features, zero-mean features
  # Cxx and Cyy are covariances from cov(Xc) and cov(Yc)
  # dz shows the desired rank of latent Z

  # NOTE: here the dimensions of Xc and Yc do not need to match
  # Note: in previous solve.w the input data was features x samples

  # Traditional CCA solution (modified from cancor function):
  nr <- nrow(Xc)
  qx <- qr(Xc)
  qy <- qr(Yc)
  dx <- qx$rank
  dy <- qy$rank

  if (dx < ncol(Xc) || dy < ncol(Yc)) {
    stop("Unable to calculate pCCA; sample covariance matrix not invertible.")
  }
  
  z <- svd(qr.qty(qx, qr.qy(qy, diag(1, nr, dy)))[1L:dx, , drop = FALSE], dx, dy)

  xcoef <- backsolve((qx$qr)[1L:dx, 1L:dx, drop = FALSE], z$u)
  ycoef <- backsolve((qy$qr)[1L:dy, 1L:dy, drop = FALSE], z$v)

  # Solve W using Archambeau06 equations
  
  # Note: only requirement for Q is that Qx%*%t(Qy) = canonical correlations
  # Q corresponds to M (dz x dz) in Bach-Jordan 2005, p.8 (before sec 4.1)
  Qx <- diag(z$d[1:dz], dz, dz) # dz x dz matrix
  #Qy <- diag(1, nrow(Qx)) # also a dz x dz matrix: identity matrix -> omit

  # ML estimates for the prob. model W:
  dz <- ifelse(is.null(dz), length(z$d), dz)
  Wx <- Cxx%*%xcoef[,1:dz]%*%Qx
  Wy <- Cyy%*%ycoef[,1:dz]#%*%Qy # Qy is identity matrix -> omit

  list(X = Wx, Y = Wy)
}

solve.archambeau <- function (X, Y, Wx, Wy, btb.x, btb.y) {

  # (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen          
  # All rights reserved. 
  # FreeBSD License (keep this notice)     

  # Use the trick introduced in Archambeau et al., ICML 2006. Robust
  # probabilistic projections, and the later correction appendix.

  # This is used to retrieve the original principal components from
  # the probabilistic CCA model solution i.e. get round the rotational
  # invariance problem

  # NOTE: Archambeau uses B in a different meaning than Bach-Jordan.                                  
  # Here denote his B with B.arch  
	
  # btb = B%*%t(B)   

  Nd <- nrow(Wx)

  Bx.arch <- Wx%*%solve(btb.x)%*%t(Wx) + diag(Nd)
  By.arch <- Wy%*%solve(btb.y)%*%t(Wy) + diag(Nd)
				
  # R is rotation matrix
  R <- eigen((diag(Nd) - solve(Bx.arch))%*%(diag(Nd) - solve(By.arch)))$vector
  Ux <- solve(cov(t(X)))%*%Wx%*%solve(matrix.sqrt((diag(Nd)-solve(Bx.arch))))%*%R
  Uy <- solve(cov(t(Y)))%*%Wy%*%solve(matrix.sqrt((diag(Nd)-solve(By.arch))))%*%R

  list(X = Ux, Y = Uy)

}
