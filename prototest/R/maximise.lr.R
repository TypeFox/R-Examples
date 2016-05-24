maximise.lr <-
function(y, U, mu=NULL, sigma=NULL, theta=NULL, exact=TRUE, verbose=FALSE, tol=10^-8, maxit=10000){
  ### first some precompuations
  n = length(y)
  M = ncol(U)
  Uy = t(U)%*%y
  y1 = sum(y)
  yPy = sum((Uy)^2)
  Py = U%*%Uy
  yP1 = sum(Py)
  
  ## for checks if whether we have specified anty parameters
  orig.mu = mu
  orig.sigma = sigma
  orig.theta = theta
  
  ### now iterate updates until convergence
  theta = 0 # have to start somewhere
  iter = 0
  while (TRUE){
    
    mu.change = 0
    sigma.change = 0
    theta.change = 0
    if (is.null(orig.mu)){ # update mu
      new.mu = update.mu (theta=theta, y1=y1, yP1=yP1, n=n)
      mu.change = abs(mu - new.mu)
      mu = new.mu
    }
    if (is.null(orig.sigma)){ # update sigma
      new.sigma = update.sigma (theta=theta, mu=mu, y=y, Py=Py)
      sigma.change = abs(sigma-new.sigma)
      sigma=new.sigma
    }
    if (is.null(orig.theta)){ # update the theta (different updates for the exact and approximate likelihoods)
      if(exact){new.theta = update.theta(mu=mu, sigma=sigma, yPy=yPy, yP1=yP1, M=M)}
      else{new.theta = update.theta.approx(mu=mu, sigma=sigma, yPy=yPy, yP1=yP1, M=M)}
      theta.change = abs(new.theta-theta)
      theta = new.theta
    }
  
    if (verbose){
      if (exact){
        ll=compute.exact.lr (theta=theta, mu=mu, sigma=sigma, y=y, Py=Py, M=M)
      }else{
        ll=compute.approx.lr (theta=theta, mu=mu, sigma=sigma, y=y, Py=Py, M=M)
      }
      print(paste("loglik: ", ll), quote=FALSE)
      print(paste("    Theta: ", theta, "; change: ", theta.change, sep=""), quote=FALSE)
      print(paste("    Mu: ", mu, "; change: ", mu.change, sep=""), quote=FALSE)
      print(paste("    Sigma: ", sigma, "; change: ", sigma.change, sep=""), quote=FALSE)
    }
    
    ## check for convergence
    if (max(mu.change, sigma.change, theta.change) < tol){break}
    
    ## count the iterations
    iter = iter + 1
    if (iter >= maxit){
      warning("exact LR stat compuation: maximum iterations reached")
      break
    }
  }
  
  if (exact){
    ll=compute.exact.lr (theta=theta, mu=mu, sigma=sigma, y=y, Py=Py, M=M)
  }else{
    ll=compute.approx.lr (theta=theta, mu=mu, sigma=sigma, y=y, Py=Py, M=M)
  }
  
  list (ll=ll, theta=theta, mu=mu, sigma=sigma)
}
