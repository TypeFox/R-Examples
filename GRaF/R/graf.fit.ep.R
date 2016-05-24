graf.fit.ep <-
  function(y, x, mn, l, wt, e = NULL, tol = 1e-6, itmax = 50, itmin = 2,
           verbose = FALSE) {
    # fit a GRaF model using expectation-propagation
    # as implemented in the gpml matlab library
    # If parallel  = TRUE, the EP uses the parallel update
    # as implemented in GPStuff
    
    # whether to use parallel EP (rather than sequential)
    parallel = FALSE
    
    if (is.vector(x)) { 
      x <- as.matrix(x)
    }
    
    # mn to probability scale
    mn <- qnorm(mn)
    n <- length(y)
    # covariance matrix
    K <- cov.SE(x1 = x, e1 = e, e2 = NULL, l = l)
    # identity matrix
    eye <- diag(n)
    
    # convert observations to +1, -1, save 0, 1 version
    oldy <- y
    y <- y * 2 - 1
    
    # initialise
    
    # \tilde{\mathbf{\nu}} = \mathbf{0}
    tnu <- rep(0, n)
    # \tilde{\mathbf{\tau}} = \mathbf{0}
    ttau <- rep(0, n)
    # \mathbf{\mu} = \mathbf{0}
    mu <- rep(0, n)
    # \Sigma = \mathbf{K}  (only used in sequential EP)
    Sigma <- K
    
    # calculate marginal negative log likelihood at ttau = tnu = mu = 0s
    z <- mu / sqrt(1 + diag(K))
    mnll <- -sum(pnorm(y * z), log.p = TRUE)

    # ~~~~~~~~~~~~~~~~~~~~~
    # set up for EP algorithm
    # set up damping factor (for parallel EP) as in GPStuff
    df <- 0.8
    
    # set old mnll to Inf & start iteration counter
    mnll_old <- Inf
    logM0_old <- logM0 <- rep(0, n)
    it <- 1
    converged <- FALSE
    
    while (!converged) {
      
      if (parallel) {
                        
        # calculate key parameters
        dSigma <- diag(Sigma)
        tau <- 1 / dSigma - ttau
        nu <- 1 / dSigma * mn - tnu
        mu <- nu / tau
        sigma2 <- 1 / tau
        
        # get marginal moments of posterior
        lis <- ep.moments(y, sigma2, mu)
        
        # recalculate parameters from these moments
        # \delta\tilde{\tau} (change in \tilde{\tau})
        dttau <- 1 / lis$sigma2hat - tau - ttau
        
        # \tilde{\tau}
        ttau <- ttau + df * dttau
        
        # \delta\tilde{\nu} (change in \tilde{\nu})
        dtnu <- 1 / lis$sigma2hat * lis$muhat - nu - tnu
        
        # \tilde{\nu}
        tnu <- tnu + df * dtnu
        
      } else {
        # otherwise use sequential EP
        
        lis <- list(Sigma = Sigma,
                    ttau = ttau,
                    tnu = tnu,
                    mu = mu)
        
        # cycle through in random order
        for (i in sample(1:n)) {
          lis <- update.ep(i, y, mn, lis)    
        } # end permuted for loop
        
        Sigma <- lis$Sigma
        ttau <- lis$ttau
        tnu <- lis$tnu
        mu <- lis$mu
        
        sigma2 <- diag(Sigma)
        
      }  # end parallel / sequential
      
      # recompute the approximate posterior parameters \Sigma and \mathbf{\mu}
      # using eq. 3.53 and eq. 3.68
      
      # sW = \tilde{S}^{\frac{1}{2}} = \sqrt{\mathbf{\tilde{\tau}}}
      sW <- sqrt(ttau)
      
      # L = cholesky(I_n  + \tilde{S}^{\frac{1}{2}} * K * \tilde{S}^{\frac{1}{2}})
      L <- chol(sW %*% t(sW) * K + eye)
      
      # V = L^T \\ \tilde{S}^{\frac{1}{2}} * K
      sWmat <- matrix(rep(sW, n), n) # byrow = TRUE?
      V <- backsolve(L, sWmat * K, transpose = TRUE)
      
      # \Sigma = \mathbf{K} - \mathbf{V}^T \mathbf{V}
      Sigma <- K - t(V) %*% V
      
      # \mathbf{\mu} = \Sigma \tilde{\mathbf{\nu}}
      mu <- Sigma %*% tnu  # + mn
      
      # calculate new mnll and assess convergence
      # compute logZ_{EP} using eq. 3.65, 3.73 and 3.74 and the existing L
      # \mathbf{\sigma}^2 = diag(\Sigma)
      sigma2 <- diag(Sigma)
      tau_n <- 1 / sigma2 - ttau
      nu_n <- mu / sigma2 - tnu + mn * tau_n
      
      z <- nu_n / tau_n / (sqrt(1 + 1 / tau_n))
      lZ <- pnorm(y * z, log.p = TRUE)
      
      # split the final equation up into 5 bits...
      mnll.a <- sum(log(diag(L))) - sum(lZ) - t(tnu) %*% Sigma %*% tnu / 2
      mnll.b <- t(nu_n - mn * tau_n)
      mnll.c <- ((ttau / tau_n * (nu_n - mn * tau_n) - 2 * tnu) / (ttau + tau_n)) / 2
      mnll.d <- sum(tnu ^ 2 / (tau_n + ttau)) / 2
      mnll.e <- sum(log(1 + ttau / tau_n)) / 2
      
      mnll <- as.numeric(mnll.a - mnll.b %*% mnll.c + mnll.d - mnll.e)
      
      # improvement in negative log marginal likelihood
      dmnll <- abs(mnll - mnll_old)
      
      # improvement in log of the 0th moment
      dlogM0 <- max(abs(logM0 - logM0_old))
      
      # both under tolerance?
      sub_tol <- dmnll < tol & dlogM0 < tol
      
      if ((sub_tol & it >= itmin) | it >= itmax) {
        # stop if there was little improvement and there have been at least
        # itmin iterations or if there have been itmax or more
        # iterations
        converged <- TRUE
      } else {
        it <- it + 1
        mnll_old <- mnll
      }
      
    } # end while loop
    
    # throw an error if the iterations maxed out before convergence
    if (it >= itmax) {
      stop(paste0('maximum number of iterations (', itmax, ') reached
                  without convergence'))
    }
    
    # calculate posterior parameters
    sW <- sqrt(ttau)
    alpha = tnu - sW * backsolve(L, backsolve(L, sW * (K %*% tnu), transpose = TRUE))
    f = crossprod(K, alpha) + mn
    
    # return relevant parameters
    return (list(y = oldy,
                 x = x,
                 MAP = f,
                 ls = l,
                 a = alpha,
                 W = ttau,
                 L = L,
                 K = K,
                 e = e,
                 obsx = x,
                 obsy = oldy,
                 mnll = mnll,
                 wt = wt))
}