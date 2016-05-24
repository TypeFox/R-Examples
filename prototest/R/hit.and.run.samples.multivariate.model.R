#### generates hit-and-run samples for the multivariate model
#### y = mu + sum_{k != test.group}theta_k.P_k(y-mu) + eps
#### conditions on the projection of y onto the column space of the nuisance (i.e. non-test.group) columns
#### distinguishes between cases where we know mu and where we do not (conditioning on additional information when mu is NULL)
#### assumes we know sigma
hit.and.run.samples.multivariate.model <-
function(x, y, groups, test.group, A, b, mu, sigma, hr.iter, hr.burn.in){
  n = length(y)
  
  ### find P2.tilde -- the projection matrix onto the nuisance columns
  nuisance.selector = groups != test.group
  if (sum(nuisance.selector) == 0){ # no columns selected in the nuisance groups
    if (is.null(mu)){
      P2 = matrix (1/n, nrow=n, ncol=n)
    }else{
      P2 = matrix (0, nrow=n, ncol=n)
    }
  }else{
    if (is.null(mu)){
      X2 = cbind(rep(1, n), x[,nuisance.selector,drop=FALSE])
    }else{
      X2 = x[,nuisance.selector,drop=FALSE]
    }
    P2 = X2%*%ginv(X2)
  }

  ### constraint matrices
  A.row.sum = apply (A, 1, sum)
  A.tilde = sigma*A%*%(diag(n) - P2)
  if (is.null(mu)){
    delta = P2%*%y
    b.tilde = b - A%*%delta
    y.init = (y - delta)/sigma
  }else{
    delta = P2%*%(y-mu)
    b.tilde = b - mu*A.row.sum - A%*%delta
    y.init = (y - mu - delta)/sigma
  }
  
  ### generate hit-and-run samples for the random part
  eps.hr = sigma*rcpp_generate_hit_and_run_samples (num_samples=hr.iter, burn_in=hr.burn.in, init_y = y.init, A=A.tilde, b=b.tilde)
  y.hr = eps.hr - P2%*%eps.hr # adjust the covariance matrix
  
  ### add back the conditioned on part
  y.hr = apply (y.hr, 2, function(col){col + delta}) # add back the projection onto the nuisance columns
  if (!is.null(mu)){y.hr = y.hr + mu} # if we have a mean, add it back
  
  return (y.hr)
}
