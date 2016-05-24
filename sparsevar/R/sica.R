# % sica.m
# %
# % It computes a minimizer of the SICA-penalized least-squares problem:
# %
# % min_beta  (2 n)^(-1) ||y - X beta||_2^2 + lambda ||rho_a(|beta|)||_1,
# %
# % where y is an n-vector of response, X is an n x p design matrix with
# % each column vector rescaled to have L2-norm n^{1/2}, lambda >= 0 is the
# % regularization parameter, and rho_a(t) = (a + 1)*t/(a + t), t >= 0, is
# % the smooth integration of counting and absolute deviation (SICA) penalty
# % (Lv and Fan, 2009) with shape parameter 0 <= a <= Infinity.
# %
# % SICA provides a family of concave penalty functions connecting the L0-
# % and L1-penalties. The L0-penalty is the target penalty function for
# % sparse recovery in linear equations, and the L1-penalty is used in
# % L1-regularization methods such as the Lasso.
# %
# % This code, which uses multi-scale stabilization, implements the SICA 
# % regularization method with the coordinate descent algorithm, following the 
# % idea of the iterative coordinate ascent (ICA) algorithm (Fan and Lv, 2011).
# %
# % References:
# %
# % 1. Fan, J. and Lv, J. (2011). Nonconcave penalized likelihood with
# % NP-dimensionality. IEEE Transactions on Information Theory 57, 5467-5484.
# %
# % 2. Lv, J. and Fan, Y. (2009). A unified approach to model selection and
# % sparse recovery using regularized least squares. The Annals of Statistics
# % 37, 3498-3528.
# %
# % Written by: Yingying Fan and Jinchi Lv, University of Southern California
# % Email: fanyingy@marshall.usc.edu
# %        jinchilv@marshall.usc.edu
# % Website: http://www-bcf.usc.edu/~fanyingy
# %          http://www-bcf.usc.edu/~jinchilv
# %
# % This version: August 1, 2011

sica <- function(X, y, a = 1e-3, lambda = 1e-2, inival = integer(), maxiter = 50, tol = 1e-4) {
  
  n <- nrow(X)
  p <- ncol(X)
  
  a <- max(1e-3, a)
  
  # % varset is an initial set of variables believed to be important and will
  # % be inlcuded for computation in each iteration
  if (length(inival) < p) {
    varset <- inival
    inival <- rep(0, p)
  } else {
    varset <- which(inival != 0) # TODO: CHECK
  }
  
  # % rescale X to make each column vector have L2-norm n^{1/2}
  Xsca <- sqrt(colSums(X^2))/sqrt(n)                         # TODO: CHECK
  X <- X / (matrix(rep(Xsca, n), n, p, byrow = TRUE))        # TODO: CHECK
  
  XXmat <- (1/n) * t(X) %*% X
  cvec <- (1/n) * t(X) %*% y
  a0 <- a
  
  # % Multi-scale stabilization
  # 
  # % First stabilization using an intermediate SICA penalty with a large shape 
  # % parameter a = 1, making the maximum concavity lambda*2*(1/a + 1/a^2) of
  # % the penalty lambda*rho_a at a low level
  if ((a0 < 1) & (4*lambda > 1e-2)) {
    a <- 1
    
    # % first round of iteration
    beta <- inival
    iter <- 1
    update <- 1
    ind <- 1:p
    
    while ((iter <= maxiter) & (update > tol)) {
      
      iter <- iter + 1
      betaold <- beta
      
      for (k in 1:length(ind)) {
        setr <- ind
        I <- setr[k]
        setr[k] <- numeric()
      }
      
      # % solve
      # %
      # % min_beta  2^(-1) (beta - z)^2 + lambda rho_a(|beta|)
      # %
      # % for scalar beta
      if (length(setr) == 0) {
        z <- cvec[I]
      } else{
        z <- (cvec[I] - XXmat[I, setr]%*%beta[setr])
      }
      beta[I] <- usica(z, a, lambda)
    }
    
    # TODO: check
    ind <- which(beta != 0)
    update <- sqrt(sum((beta - betaold)^2))
    setr <- setdiff(1:p, ind)
    
    if (length(setr) == 0) {
      resc <- 0
    } else {
      resc <- abs(cvec(setr) - XXmat(setr, ind)*beta(ind))
    }
    
    indm <- which(resc > lambda*(1 + a^(-1)))
    ind <- union(c(setr(indm), ind), varset)    
  }
  
  # TODO: check this two lines
  varset <- union(which(beta != 0), varset)
  inival <- beta
  
  # % Second stabilization using an intermediate SICA penalty with a relatively 
  # % large shape parameter a = 1/3, making the maximum concavity lambda*2*(1/a + 1/a^2) 
  # % of the penalty lambda*rho_a at a relatively low level
  
  if ((a0 < 1/3) & (lambda*24 > 1e-2)) {
    
    a <- 1/3
    
    # % first round of iteration
    beta <- inival
    iter <- 1
    update <- 1
    ind <- 1:p
    
    while ((iter <= maxiter) & (update > tol)) { 
      
      iter <- iter + 1
      betaold <- beta
      
      for (k in 1:length(ind)) {
        setr <- ind
        I <- setr[k]
        setr[k] <- numeric(0)
        
        # % solve
        # %
        # % min_beta  2^(-1) (beta - z)^2 + lambda rho_a(|beta|)
        # %
        # % for scalar beta
        
        if (length(setr) == 0) {
          z <- cvec[I]
        } else {
          z <- (cvec[I] - XXmat(I, setr)*beta(setr))
        }
        
        beta[I] <- usica(z, a, lambda)
      }
      
      ind <- which(beta != 0)
      # TODO: check substitute for norm
      # update <- norm(beta - betaold);
      update <- sqrt(sum((beta - betaold)^2))
      
      setr <- setdiff(1:p, ind)
      
      if (length(setr) == 0){
        resc <- 0
      } else {
        resc <- abs(cvec(setr) - XXmat(setr, ind)*beta(ind))
      }
      
      indm <- which(resc > lambda*(1 + a^(-1)))
      ind <- union(c(t(setr(indm)), ind), varset)
    }
    
    varset <- union(which(beta != 0), varset)
    inival <- beta
  }
  
  ########################################################################
  # CHECK & CONVERT
  ########################################################################
  
  ## % Third stabilization using an intermediate SICA penalty with a relatively 
  ## % large shape parameter a = 0.1, making the maximum concavity lambda*2*(1/a + 1/a^2) 
  ## % of the penalty lambda*rho_a at a relatively low level
  
  if ((a0 < 0.1) & (lambda*220 > 1e-2)) {
    a <- 0.1
    
    # % first round of iteration
    beta <- inival
    iter <- 1
    update <- 1
    ind <- 1:p
    
    while ((iter <= maxiter) & (update > tol)) {
      iter <- iter + 1
      betaold <- beta
      
      for (k in 1:length(ind)) {
        setr <- ind
        I <- setr[k]
        setr[k] <- numeric(0)
        
        ## % solve
        ## %
        ## % min_beta  2^(-1) (beta - z)^2 + lambda rho_a(|beta|)
        ## %
        ## % for scalar beta
        if (length(setr) > 0) {
          z <- cvec[I]
        } else {
          z <- (cvec[I] - XXmat(I, setr)*beta[setr]);
        }
        beta[I] <- usica(z, a, lambda)
      }
      
      ind <- which(beta != 0)
      update <- sqrt(sum((beta - betaold)^2))
      
      setr <- setdiff(1:p, ind)
      if (length(setr) == 0) {
        resc <- 0
      } else {
        resc  <- abs(cvec(setr) - XXmat(setr, ind)*beta(ind))
      }
      
      indm <- which(resc > lambda*(1 + a^(-1)))
      ind <- union(c(t(setr(indm)), ind), varset)
    }
    
    varset <- union(which(beta != 0), varset)
    inival <- beta
  }
  
  # % Final solution
  a <- a0
  # % first round of iteration
  beta <- inival
  iter <- 1
  update <- 1
  ind <- 1:p
  
  while ((iter <= maxiter) & (update > tol)) {
    
    iter <- iter + 1
    betaold <- beta
    
    for (k in 1:length(ind)) {
      setr <- ind
      I <- setr[k]
      setr[k] <- numeric(0)
      
      ## % solve
      ## %
      ## % min_beta  2^(-1) (beta - z)^2 + lambda rho_a(|beta|)
      ## %
      ## % for scalar beta
      if (length(setr) == 0) {
        z <- cvec[I]
      } else {            
        z <- (cvec[I] - XXmat(I, setr)*beta(setr))
      }
      
      beta[I] <- usica(z, a, lambda)
    }
    
    ind <- utils::find(beta)
    update <- norm(beta - betaold)
    
    setr <- setdiff(1:p, ind)
    
    if (length(setr) == 0) {
      resc <- 0
    } else {
      resc <- abs(cvec(setr) - XXmat(setr, ind)*beta(ind))
    }
    
    indm <- utils::find(resc > lambda*(1 + a^(-1)))
    ind <- union(c(t(setr(indm)), ind), varset)
  }
  
  # % rescale beta vector to original scale
  beta <- beta/t(Xsca)
  
  # % If the size of selected model exceeds n/2, display a warning message
  if (sum(beta != 0) > n/2) {
    cat(" ")
    cat("Warning: The solution found is nonsparse and may be inaccurate. Try a larger lambda!")
    cat(" ")
  }
  
}

# % usica.m
# %
# % It computes the global minimizer of the univariate SICA-penalized
# % least-squares problem:
# %
# % min_beta  (2 Lam)^(-1) (beta - z)^2 + lambda rho_a(|beta|),
# %
# % where Lam > 0 (default value 1 in linear models), z is a scalar, lambda
# % >= 0 is the regularization parameter, and rho_a(t) = (a + 1)*t/(a + t),
# % t >= 0, is the smooth integration of counting and absolute deviation
# % (SICA) penalty (Lv and Fan, 2009) with shape parameter 0 <= a <=
# % Infinity. This minimization problem involves solving a cubic equation,
# % and thus the SICA solution in one dimension admits an analytical form.
# %
# % SICA provides a family of concave penalty functions connecting the L0-
# % and L1-penalties. The L0-penalty is the target penalty function for
# % sparse recovery in linear equations, and the L1-penalty is used in
# % L1-regularization methods such as the Lasso.
# %
# % This code is a key step in the coordinate descent algorithm for 
# % implementing the SICA regularization method, following the idea of the 
# % iterative coordinate ascent (ICA) algorithm (Fan and Lv, 2011).
# %
# % References:
# %
# % 1. Fan, J. and Lv, J. (2011). Nonconcave penalized likelihood with
# % NP-dimensionality. IEEE Transactions on Information Theory 57, 5467-5484.
# %
# % 2. Lv, J. and Fan, Y. (2009). A unified approach to model selection and
# % sparse recovery using regularized least squares. The Annals of Statistics
# % 37, 3498-3528.
# %
# % Written by: Yingying Fan and Jinchi Lv, University of Southern California
# % Email: fanyingy@marshall.usc.edu
# %        jinchilv@marshall.usc.edu
# % Website: http://www-bcf.usc.edu/~fanyingy
# %          http://www-bcf.usc.edu/~jinchilv
# %
# % This version: August 1, 2011
# %

######################################################
# TESTED AGAINST MATLAB OUTPUTS: SEEMS TO BE WORKING #
######################################################

usica <- function(z, a, lambda, Lam = 1) {
  
  #sicap = @(t, a) (a + 1)*abs(t)./(a + abs(t))
  a <- max(1e-3, a)
  
  if (a > 1e3){
    beta <- sign(z)*max(0, abs(z) - Lam*lambda)
  } else {
    b <- 2*a - abs(z)
    c <- a^2 - 2*a*abs(z)
    d <- Lam*lambda*a*(a + 1) - a^2*abs(z)
    q <- (3*c - b^2)/9
    r <- (9*b*c - 27*d - 2*b^3)/54
    D <- q^3 + r^2
    
    if (D > 0){
      beta <- 0
    } else { 
      tvec <- sort(Re(polyroot(c(d, c, b, 1))));
      z0 <- Lam*lambda*(1 + a^(-1));
      if (abs(z) >= z0){
        beta <- sign(z)*tvec[3];
      } else if ((tvec[2] >= 0) & (tvec[3] <= abs(z)) & (tvec[2] < tvec[3]) & ((2*Lam)^(-1)*z^2 > (2*Lam)^(-1)*(abs(z) - tvec[3])^2 + lambda*sicap(tvec[3], a) + 1e-8)) {
        beta <- sign(z)*tvec[3]
      } else {
        beta <- 0
      }
    }
    
  }
  return(beta)
}

sicap <- function(t, a) {
  
  return((a + 1)*abs(t)/(a + abs(t)))
  
}
