fitfrail.fit <- function(x, y, cluster, init.beta, init.theta, frailty, 
                         control, rownames, weights) {
  
  # TODO: error check for number of frailty distr params
  if (missing(init.theta) || length(init.theta) != length(init.frailty[[frailty]]))
    stop("Wrong length for frailty distribution inital values")
    
  if (!missing(init.beta) && length(init.beta) != ncol(x))
    stop("Wrong length for coefficient inital values")
  
  # VARS acts as a cache and holds all variables needed for computations
  # Some variables are shared across functions, e.g. score and jacobian fns
  # Modifying the execution environment of fitfrail.fit make the function much 
  # more complex, but avoids repeating expensive computations
  VARS <- environment()
  
  # Size of beta, theta, and gamma = c(beta, theta)
  n.beta <- length(init.beta)
  n.theta <- length(init.theta)
  n.gamma <- n.beta + n.theta
  init.gamma <- c(init.beta, init.theta)
  
  # time and failure indicator
  time <- y[,1]
  status <- y[,2]
  
  time_sorted_idx <- order(time)
  time_steps <- c(0, time[time_sorted_idx])
  # time_steps <- c(0, sort(unique(time)))
  k_tau <- length(time_steps) # index of the last time step
  
  # d_[k] is the number of failures at time t_k
#   d_ <- vapply(seq_along(time_steps), function(k) {
#     sum(status[time == time_steps[k]])
#   }, 0)
  d_ <- c(0, unname(status[time_sorted_idx]))
  
  # Variables indexed by [[i]][j,] where i is cluster name and j is member
  X_ <- split.data.frame(x, cluster) # Covariates
  I_ <- split(status, cluster) # Failure indicator
  
  # Helper fn to get the rank, where duplicates have the same rank
  increasing_rank <- function(d) {
    j <- unique(sort(d));
    return(sapply(d,function(dd) which(dd==j)));
  }
  
  # Index of the observed time, K_[[i]][j]
  # K_ <- split(1 + increasing_rank(time), cluster) 
  K_ <- split(1 + rank(time, ties.method="first"), cluster)
  
  # Cluster names and sizes
  cluster_names <- levels(cluster)
  cluster_sizes <- table(cluster)
  names(cluster_names) <- cluster_names
  n.clusters <- length(cluster_names)
  
  if (is.null(weights)) {
    weights <- rep(1, n.clusters)
  } else {
    stopifnot(length(weights) == n.clusters)
  }
  
  # Convenience functions to apply over the clusters and members
  clustapply <- function(fn, FUN.VALUE) {
    lapply(cluster_names, function(i) {
      t(vapply(1:cluster_sizes[[i]], function(j) {
        fn(i, j)
      }, FUN.VALUE=FUN.VALUE))})
  }
  
  clustlapply <- function(fn) {
    lapply(cluster_names, function(i) {
      lapply(1:cluster_sizes[[i]], function(j) {
        fn(i, j)
      })})
  }
  
  # N_[[i]][j,k] == 1 if ij failed on/before time t_k
  N_ <- clustapply(function(i, j) {
    as.numeric((I_[[i]][j] > 0) & (K_[[i]][j] <= 1:k_tau))
    }, rep(0, k_tau))
  
  # N_dot_[[i]][k] is the number of individuals in cluster i that failed at time <= t_k
  N_dot_ <- lapply(N_, colSums)
  
  # N_tilde_[[i]][j,k] == 1 if ij was observed at time <= t_k
  N_tilde_ <- clustapply(function(i, j) {
    as.numeric(K_[[i]][j] <= 1:k_tau)
    }, rep(0, k_tau))
  
  # Y_[[i]][j,k] == 1 if ij failed at time >= t_k
  Y_ <- clustapply(function(i, j) {
      as.numeric(K_[[i]][j] >= 1:k_tau)
    }, rep(0, k_tau))
  
  iter <- 0
  trace <- matrix(nrow=0, ncol=n.gamma+2)
  colnames(trace) <- c("Iteration",
                       paste("beta.", names(init.beta), sep=""),
                       paste("theta.", 1:n.theta, sep=""), 
                       "loglik")
  
  # The objective function for either score equations or loglikelihood optimization
  # gamma is c(beta, theta), several global variables are updated since these are
  # reused in other places (jacobian, covar matrix)
  fit_fn <- function(hat.gamma) {
    # Update the current estimate. Everything below depends on this
    VARS$hat.gamma <- hat.gamma
    
    # Modify the execution environment of fitfrail.fit
    with(VARS, {
      if (iter >= control$maxit) {
        if (control$fitmethod == "loglik") {
          return(loglik)
        } else if (control$fitmethod == "score") {
          return(U_)
        }
      }
      iter <- iter + 1
      
      hat.beta <- hat.gamma[1:n.beta]
      hat.theta <- hat.gamma[(n.beta+1):(n.beta+n.theta)]
      
      R_ <- clustapply(function(i, j) {
        exp(sum(hat.beta * X_[[i]][j,])) * Y_[[i]][j,]
      }, rep(0, k_tau))
      R_dot_ <- lapply(R_, colSums)
      
      R_star <- clustapply(function(i, j) {
        exp(sum(hat.beta * X_[[i]][j,]))
      }, 0)
      
      bh_ <- bh(d_, R_star, K_, Y_, N_, N_dot_, hat.beta, hat.theta, frailty, weights,
                control$int.abstol, control$int.reltol, control$int.maxit)
      
      H_ <- bh_$H_
      H_dot_ <- bh_$H_dot
      lambda <- bh_$lambda
      Lambda <- bh_$Lambda
      psi_ <- bh_$psi_
      phi_1_ <- bh_$phi_1_
      phi_2_ <- bh_$phi_2_
      
      phi_prime_1_ <- lapply(1:n.theta, function(theta_idx) 
        phi_prime_k(1, theta_idx, N_dot_, H_dot_, hat.theta, frailty,
                    control$int.abstol, control$int.reltol, control$int.maxit))
      
      loglik_vec <- loglikelihood(X_, K_, I_, phi_1_, lambda, hat.beta)
      loglik <- sum(weights * loglik_vec)
      
      if (control$verbose) {
        if (iter == 1) {
          cat(c("Iter", 
                paste("beta.", names(init.beta), sep=""),
                paste("theta.", 1:n.theta, sep=""), 
                "loglik",
                "\n"), sep="\t")
        }
        cat(c(iter, format(round(c(hat.beta, hat.theta, loglik), 
                           4), nsmall=4, trim=TRUE), "\n"), sep="\t")
      }
      
      xi_beta_ <- matrix(vapply(1:n.beta, function(r) {
        xi_beta(X_, I_, H_, psi_, r)
      }, rep(0, n.clusters)), nrow=n.clusters, ncol=n.beta)
      
      xi_theta_ <- matrix(vapply(1:n.theta, function(r) {
        xi_theta(phi_1_, phi_prime_1_[[r]], r)
      }, rep(0, n.clusters)), nrow=n.clusters, ncol=n.theta)
      
      # xi_[i,r], where i is cluster, r is param idx
      # Update globally, this is used in the jacobian and covariance matrix
      xi_ <- cbind(xi_beta_, xi_theta_)
      U_ <- colMeans(weights * xi_)
      
      trace <- rbind(trace, c(iter, unname(hat.gamma), loglik))
      
      if (control$fitmethod == "loglik") {
        loglik
      } else if (control$fitmethod == "score") {
        U_
      }
  })}
  
  # Already computed
  loglik_jacobian <- function(gamma=NULL) with(VARS, {
    U_
  })
  
  score_jacobian <- function(gamma=NULL) with(VARS, {
    phi_3_ <- phi_k(3, N_dot_, H_dot_, hat.theta, frailty,
                    control$int.abstol, control$int.reltol, control$int.maxit)
    
    phi_prime_2_ <- lapply(1:n.theta, function(theta_idx) 
      phi_prime_k(2, theta_idx, N_dot_, H_dot_, hat.theta, frailty,
                  control$int.abstol, control$int.reltol, control$int.maxit))
    
    phi_prime_prime_1_ <- lapply(1:n.theta, function(theta_idx_1) {
      lapply(1:n.theta, function(theta_idx_2) {
        phi_prime_prime_k(1, theta_idx_1, theta_idx_2, 
                          N_dot_, H_dot_, hat.theta, frailty, k_tau,
                          control$int.abstol, control$int.reltol, control$int.maxit)})})
    
    dH_dbeta_ <- lapply(1:n.beta, function(s) {
      dH_dbeta(s, d_, X_, K_, R_, R_dot_, R_star,
               phi_1_, phi_2_, phi_3_, 
               Lambda, lambda,
               hat.beta, hat.theta, frailty)
    })
    
    dH_dtheta_ <- lapply(1:n.theta, function(s) {
      dH_dtheta(d_, X_, K_, R_, R_dot_, R_star,
                phi_1_, phi_2_, phi_3_,
                phi_prime_1_[[s]], phi_prime_2_[[s]],
                Lambda, lambda, hat.beta)
    })
    
    jacobian_ls <- function(l, s) {
      if (l <= n.beta && s <= n.beta) {
        jacobian_beta_beta(l, 
                           X_, K_, H_, 
                           phi_1_, phi_2_, phi_3_,
                           dH_dbeta_[[s]]$dH_dbeta_, 
                           dH_dbeta_[[s]]$dH_dot_dbeta_)
      } else if (l <= n.beta && s > n.beta) {
        theta_idx <- s - n.beta
        jacobian_beta_theta(l,
                            X_, K_,  H_, 
                            phi_1_, phi_2_, phi_3_,
                            phi_prime_1_[[theta_idx]], phi_prime_2_[[theta_idx]],
                            dH_dtheta_[[theta_idx]]$dH_dtheta_, 
                            dH_dtheta_[[theta_idx]]$dH_dot_dtheta_)
      } else if (l > n.beta && s <= n.beta) {
        theta_idx <- l - n.beta
        jacobian_theta_beta(phi_1_, phi_2_,
                            phi_prime_1_[[theta_idx]], phi_prime_2_[[theta_idx]],
                            dH_dbeta_[[s]]$dH_dot_dbeta_)
      } else if (l > n.beta && s > n.beta) {
        theta_idx_1 <- l - n.beta
        theta_idx_2 <- s - n.beta
        jacobian_theta_theta(phi_1_, phi_2_,
                             phi_prime_1_[[theta_idx_1]], phi_prime_2_[[theta_idx_1]], 
                             phi_prime_prime_1_[[theta_idx_1]][[theta_idx_2]], 
                             dH_dtheta_[[theta_idx_2]]$dH_dot_dtheta_)
      }
    }
    
    # Build the jacobian matrix from the respective components
    outer(1:n.gamma, 1:n.gamma, Vectorize(function(l, s) jacobian_ls(l, s)))
  })
  
  start.time <- Sys.time()
  # The actual optimization takes place here
  if (control$fitmethod == "loglik") {
    fitter <- optim(init.gamma, fit_fn,
                    gr=loglik_jacobian,
                    lower=c(rep(-Inf, VARS$n.beta), lb.frailty[[frailty]]), 
                    upper=c(rep(Inf, VARS$n.beta),  ub.frailty[[frailty]]), 
                    method="L-BFGS-B",
                    control=list(factr=control$reltol/.Machine$double.eps, 
                                 pgtol=0, fnscale=-1)
                    )
    hat.gamma <- fitter$par
  } else if (control$fitmethod == "score") {
    fitter <- nleqslv(init.gamma, fit_fn, 
                      control=list(xtol=control$reltol, ftol=control$abstol, btol=1e-3,
                                   allowSingular=TRUE),
                      method="Newton",
                      jac=score_jacobian,
                      jacobian=TRUE
                      )
    hat.gamma <- fitter$x
  }
  fit.time <- Sys.time() - start.time
  
  if (control$verbose)
    cat("Converged after", iter, "iterations\n")
  
  if (iter >= control$maxit) {
    warning("Maximum number of iterations reached before convergence.")
  }
  
  trace <- data.frame(trace)
  
  hat.beta <- hat.gamma[1:n.beta]
  hat.theta <- hat.gamma[(n.beta+1):(n.gamma)]
  
  # Unique time steps where failures occur
  Lambda.df.cens <- aggregate(VARS$lambda, list(time_steps), sum)
  names(Lambda.df.cens) <- c("time","lambda")
  Lambda.df.cens <- data.frame(time=Lambda.df.cens$time, Lambda=cumsum(Lambda.df.cens$lambda))
  Lambda.df <- Lambda.df.cens[duplicated(Lambda.df.cens$Lambda)==FALSE,]
  
  Lambda.fun <- Vectorize(function(t) {
    if (t <= 0) {
      return(0);
    }
    Lambda.df$Lambda[sum(t >= Lambda.df$time)]
  })
  
  list(beta=hat.beta,
       theta=setNames(hat.theta, paste("theta.", 1:length(hat.theta), sep="")),
       Lambda=Lambda.df,
       Lambda.all=Lambda.df.cens,
       Lambda.fun=Lambda.fun,
       frailty=frailty,
       frailty.variance=vfrailty[[frailty]](hat.theta),
       loglik=VARS$loglik,
       iter=iter,
       fitter=fitter,
       fitmethod=control$fitmethod,
       n.clusters=n.clusters,
       trace=trace,
       fit.time=fit.time,
       # Keep the execution environment, needed for vcov
       VARS=VARS
      )
}