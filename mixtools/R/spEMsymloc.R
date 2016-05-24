## EM-like algorithm for a nonparametric univariate mixture model with
## symmetric components from a location family
spEMsymloc <- function(x, mu0, bw = bw.nrd0(x), h=bw, eps = 1e-8, maxiter=100,
                        stochastic = FALSE, verbose = FALSE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  n <- length(x)
  if (length(mu0)>1) m <- length(mu0) # mu0=centers
  	else m <- mu0 # when mu0=number of clusters
  z.hat <- matrix(0, nrow=n, ncol=m)
  fkernel <- matrix(0, nrow=n, ncol=m)
  tt0 <- proc.time()
  lambda <- rep(1/m, m)
  kmeans <- kmeans(x, mu0)
  for(j in 1:m) {
    z.hat[kmeans$cluster==j, j] <- 1
  }
  iter <- 0
  if (stochastic) {
    sumpost <- matrix(0, n, m)
  }
  finished <- FALSE
  lambda <- mu <- matrix(0,maxiter,m)
  while (!finished) {
  #while (max(abs(change)) > eps & iter < maxiter) {
    iter <- iter + 1
    t0 <- proc.time()

    ## M-Step
    lambda[iter,] <- colMeans(z.hat) 
    mu[iter,] <- apply(sweep(z.hat, 1, x, "*"), 2, mean)/lambda[iter,]

    ## density estimation step
    if(stochastic){
      z <- t(apply(z.hat, 1, function(prob) rmultinom(1, 1, prob)))
      ans <- .C("KDEsymloc", n=as.integer(n), m=as.integer(m),
                mu=as.double(mu[iter,]), x=as.double(x), bw=as.double(bw),
                z=as.double(z), f = double(n*m),
                PACKAGE="mixtools")
    } else {
      ans <- .C("KDEsymloc", n=as.integer(n), m=as.integer(m),
                mu=as.double(mu[iter,]), x=as.double(x), bw=as.double(bw),
                z=as.double(z.hat), f = double(n*m),                     
                PACKAGE="mixtools")
    }

    fkernel <- matrix(ans$f, ncol=m)
    lambda.f <- sweep(fkernel, 2, lambda[iter,], "*")

    ## E-step (for next iteration)
    z.hat <- lambda.f/rowSums(lambda.f)    
    finished <- iter >= maxiter
    if (stochastic) {
      sumpost <- sumpost + z.hat
    } else if (iter>1) { # This convergence criterion is too simplistic:
      change <- c(lambda[iter,] - lambda[iter-1,], mu[iter,]-mu[iter-1,])
      finished <- finished | (max(abs(change)) < eps)
    }
    if (verbose) {
      t1 <- proc.time()
      cat("iteration ", iter, "  lambda ", round(lambda[iter,], 4), 
          "  mu ", round(mu[iter,], 4))
      cat(" time", (t1 - t0)[3], "\n")
    }
  }
  if (verbose) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter,], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }
  if(stochastic) {
    return(structure(list(data=x, posteriors=sumpost/iter, lambda=lambda,
                          bandwidth=bw, lambdahat=colMeans(lambda), 
                        mu = mu, muhat = colMeans(mu), symmetric=TRUE),
              class="npEM"))
  } else {
    return(structure(list(data=x, posteriors=z.hat, lambda=lambda[1:iter,],
                          bandwidth=bw, lambdahat=lambda[iter,], 
                          mu = mu[1:iter,], muhat = mu[iter,], symmetric=TRUE),
                        class="npEM"))
  } 
}

