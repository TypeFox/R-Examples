vcov.fitfrail <- function(object, boot=FALSE, B=100,
                          Lambda.times=NULL, # The time points of the CBH to use
                          cores=0, ...) {
  fit <- object
  Call <- match.call()
  
  # Use the time points where CBH increases
  if (is.null(Lambda.times)) {
    Lambda.times <- fit$Lambda$time
  }
  
  # CBH variance not supportted in the estimator yet
  if (!boot) {
    Lambda.times = NULL
  }
  
  # If V has already been computed the same way and matches the size we expect
  if (!is.null(fit$VARS[["COV"]]) && !is.null(fit$VARS[["COV.call"]])) {
    new.Call <- fit$VARS[["COV.call"]]
    cache.Call <- Call
    
    # Ignore the fit object and number of cores used
    new.Call$object <- NULL
    new.Call$cores <- NULL
    cache.Call$object <- NULL
    cache.Call$cores <- NULL
    
    if(identical(new.Call, cache.Call))
      return(fit$VARS[["COV"]])
  }
  
  # Weighted bootstrap covariance
  vcov.boot <- function() {
    fn <- function(s) {
      set.seed(s)
      
      weights <- rexp(fit$n.clusters)
      weights <- weights/mean(weights)
      
      new.call <- fit$call
      new.call$weights <- weights
      new.fit <- eval(new.call)
      
      c(new.fit$beta,
        new.fit$theta,
        setNames(new.fit$Lambda.fun(Lambda.times), 
                 paste("Lambda.", format(Lambda.times, nsmall=2), sep="")))
    }
    
    hats <- t(simplify2array(plapply(B, fn, cores=cores)))
    cov(hats)
  }
  
  # Covariance conistent estimator
  vcov.estimator <- function() with(fit$VARS, {
    
    # jacobian is computed first, sets up some vars that are needed in steps 1-3
    jac <- score_jacobian()
    
    ################################################################## Step I
    
    V <- Reduce('+', lapply(1:n.clusters, function(i) {
      xi_[i,] %*% t(xi_[i,])
    }))/n.clusters
    # V <- (t(xi_) %*% xi_)/n.clusters
    
    ################################################################## Step II
    
    # Split Q up into hat.beta and hat.theta components, then create a Q_ function
    # that accesses component by index r
    Q_beta_ <- lapply(1:n.beta, function(r) {
      Q_beta(X_, K_, H_, R_star, phi_1_, phi_2_, phi_3_, r)
    })
    
    Q_theta_ <- lapply((n.beta+1):(n.beta+n.theta), function(r) {
      Q_theta(H_, R_star, phi_1_, phi_2_,
              phi_prime_1_[[r - n.beta]], phi_prime_2_[[r - n.beta]])
    })
    
    # Q_[[r]] has the same shape as H_dot_
    Q_ <- c(Q_beta_, Q_theta_)
    
    # calligraphic Y, Ycal_[t]
    Ycal_ <- Ycal(X_, R_star, Y_, psi_, hat.beta)
    
    # eta_[t]
    eta_ <- eta(phi_1_, phi_2_, phi_3_)
    
    # Upsilon
    Upsilon_ <- Upsilon(X_, R_star, K_, R_dot_, eta_, Ycal_, hat.beta)
    
    # Omega_[[i]][[j, t]
    # Compute everything
    Omega_ <- Omega(X_, R_star, N_, R_dot_, eta_, Ycal_, hat.beta)
    
    # p_hat_[t]
    p_hat_ <- p_hat(I_, Upsilon_, Omega_, N_tilde_)
    
    # pi_[[r]][s]
    pi_ <- lapply(1:n.gamma, function(r) {
      pi_r(Q_[[r]], N_tilde_, p_hat_)
    })
    
    G <- outer(1:n.gamma, 1:n.gamma, Vectorize(function(r, l) {
      G_rl(pi_[[r]], pi_[[l]], p_hat_, Ycal_, N_)
    }))
    
    ################################################################## Step III
    
    # M_hat_[[i]][j,s]
    M_hat_ <- M_hat(X_, R_star, N_, Y_, psi_, hat.beta, Lambda)
    
    # u_star_[i,r]
    u_star_ <- u_star(pi_, p_hat_, Ycal_, M_hat_)
    
    C <- outer(1:n.gamma, 1:n.gamma, Vectorize(function(r, l) {
      sum(vapply(1:n.clusters, function(i) {
        xi_[i,r] * u_star_[i,l] + xi_[i,l] * u_star_[i,r]
      }, 0))
    }))/n.clusters
    
    ################################################################## Step IV
    D_inv <- solve(jac)
    
    ################################################################## Step V
    
    (D_inv %*% (V + G + C) %*% t(D_inv))/n.clusters
  })
  
  if (boot) {
    COV <- vcov.boot()
  } else {
    COV <- vcov.estimator()
  }
  
  param.names <- c(names(fit$beta), names(fit$theta))
  
  if (length(Lambda.times) > 0) {
    param.names <- c(param.names, paste("Lambda.", format(Lambda.times, nsmall=2), sep=""))
  }
  
  rownames(COV) <- param.names
  colnames(COV) <- param.names
  
  # Cache the value for quick access later
  fit$VARS[["COV"]] <- COV
  fit$VARS[["COV.call"]] <- Call
  
  COV
}