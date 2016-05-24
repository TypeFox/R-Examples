#### computes the exact or approx likelihood ratio statistic (user specifies which) for a given y and selected columns x
#### uses compiled cpp code to compute the maximum likelihood estimate
#### y is a MATRIX with different replications contained in the columns
compute.lr.stat.multi <-  
function(x, y, groups, test.group, mu, sigma, exact=TRUE, verbose=FALSE, tol=10^-8){
  ### get the SVD of the individual group matrices
  V = do.call(cbind, lapply (unique(groups), function(g){
    svd(x[, groups==g, drop=FALSE])$u
  }))
  
  ### initial computations
  n = nrow (y)
  M = ncol (x)
  K = length (unique(groups))
  K.min.test = length(unique(groups[groups != test.group]))
  
  ### see whether any columns in the nuisance groups are selected
  ### if not, it becomes a univariate test
  if (K.min.test == 0){# univariate test
    return (rcpp_compute_lr_stat (U=V[,groups==test.group, drop=FALSE], y=y, mu=ifelse(is.null(mu), Inf, mu), sigma=sigma, exact=exact, verbose=verbose, tol=tol, maxit=10000))
  }
  
  ### ASSUME FROM HERE THAT THERE ARE SELECTED COLUMNS IN THE NUISANCE GROUPS
  if (is.null(mu)){ # we are estimating mu
    if (exact){
      init.theta = matrix (0, nrow=K, ncol=ncol(y))
      init.theta.0 = matrix (0, nrow=K.min.test, ncol=ncol(y))
      max.ll.obj = rcpp_maximise_likelihood(init_theta=init.theta, y=y, V=V, groups=groups, mu=rep(Inf, ncol(y)), sigma=sigma, sm_inv=TRUE)
      max.ll.obj.0 = rcpp_maximise_likelihood(init_theta=init.theta.0, y=y, V=V[,groups!=test.group], groups=groups[groups!=test.group], mu=max.ll.obj[K+2,], sigma=sigma, sm_inv=TRUE) # use the estimated mu from the unconstrained model when estimting the constrained model -- forces this maximum loglik to be smaller than the previous, so we dont have those pesky negative values anymore
    }else{
      max.ll.obj = rcpp_maximise_approx_likelihood(y_mat=y, V=V, groups=groups, mu=rep(Inf, ncol(y)), sigma=sigma)
      max.ll.obj.0 = rcpp_maximise_approx_likelihood(y_mat=y, V=V[,groups!=test.group], groups=groups[groups!=test.group], mu=max.ll.obj[K+2,], sigma=sigma) # use the estimated mu from the unconstrained model when estimting the constrained model -- forces this maximum loglik to be smaller than the previous, so we dont have those pesky negative values anymore
    }
  }else{
    if (exact){
      init.theta = matrix (0, nrow=K, ncol=ncol(y))
      init.theta.0 = matrix (0, nrow=K.min.test, ncol=ncol(y))
      max.ll.obj = rcpp_maximise_likelihood(init_theta=init.theta, y=y, V=V, groups=groups, mu=rep(mu, ncol(y)), sigma=sigma, sm_inv=TRUE)
      max.ll.obj.0 = rcpp_maximise_likelihood(init_theta=init.theta.0, y=y, V=V[,groups!=test.group], groups=groups[groups!=test.group], mu=rep(mu, ncol(y)), sigma=sigma, sm_inv=TRUE)
    }else{
      max.ll.obj = rcpp_maximise_approx_likelihood(y_mat=y, V=V, groups=groups, mu=rep(mu, ncol(y)), sigma=sigma)
      max.ll.obj.0 = rcpp_maximise_approx_likelihood(y_mat=y, V=V[,groups!=test.group], groups=groups[groups!=test.group], mu=rep(mu, ncol(y)), sigma=sigma)
    }
  }
  
  return (2*(max.ll.obj[1,] - max.ll.obj.0[1,]))
}