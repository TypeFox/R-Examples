## EM-like algorithm for a location-scale mixture model with
## independent repeated measures -- some ID, some not --
## where each block has its own location and scale but otherwise
## all blocks (within a component or globally, depending) have the same
## shape.

## Correction:  For now, this algorithm only implements model (17) in 
## Benaglia et al -- in other words, each component and block has exactly
## the same shape and they differ only by location and scale.
spEM <- function(x, mu0, blockid = 1:ncol(x),
                 bw=bw.nrd0(as.vector(as.matrix(x))), constbw = TRUE,
                 h=bw, eps=1e-8,
                 maxiter=500, stochastic = FALSE, verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject
  u <- match(blockid, unique(blockid)) # convert blockid to integers 1, 2, ...
  if (is.matrix(mu0)) 
    m <- dim(mu0)[1]  # mu0=centers
  else 
    m <- mu0  # mu0=number of clusters
  z.hat <- matrix(0, nrow = n, ncol = m)
  tt0 <-  proc.time() # for total time

  ## Initial Values
  if(m == 1) z.hat <- matrix(1, nrow = n, ncol = m)
  else{
    kmeans <- kmeans(x, mu0)
    for(j in 1:m)
      z.hat[kmeans$cluster==j, j] <- 1
  }
  iter <- 0
  if (stochastic) {
    sumpost <- matrix(0, n, m)
  }
  finished <- FALSE
  lambda <- matrix(0, nrow = maxiter, ncol = m)
  mu <- sigma <- array(0, dim=c(maxiter, m, max(u)))
  stackedx <- x[rep(1:n,m),]
  loglik <- NULL
  #  Cfunction <- ifelse(constbw, "KDErepeated", "KDErepeatedbw")
  Cfunction <- "KDElocscale"
  
  while (!finished) {
    iter <- iter + 1
    bw.old <- bw
    t0 <- proc.time()
    
    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat is in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)
    
    if (stochastic) {
      z <- t(apply(z.hat, 1, function(prob) rmultinom(1, 1, prob)))
      cs <- colSums(z)
      z.tmp <- sweep(z, 2, cs, "/")
      z.tmp[, cs==0] <- 1/NROW(z.tmp) # Just in case
    }
    else {
      cs <- colSums(z.hat)
      z.tmp <- sweep(z.hat, 2, cs, "/")
      z.tmp[, cs==0] <- 1/NROW(z.tmp) # Just in case
    }
    h <- bw 
    
    ## More M-step (means and std devs are location / scale params)
    for (k in 1:max(u)) { # k is the block number
      r2 <- sum(u == k)
      x2 <- x[, u==k] # Subset of data belonging to kth block (n x r2 matrix)
      mu[iter, , k] <- as.vector(rowMeans(t(z.tmp) %*% x2))
      for (j in 1:m)
        sigma[iter, j, k] <- sqrt(sum(z.tmp[,j] * (x2-mu[iter, j, k])^2)/r2)
    }
    
    ## density estimation step
    if (!constbw) {
      wts <- rep(as.vector(z.tmp),r)
      scaledx <- as.vector((stackedx - mu[iter, rep(1:m, each=n), u])/
                           sigma[iter, rep(1:m, each=n), u])
      h <- bw <- 0.9 / (n*r)^(1/5) * min(1, wiqr<-wIQR(wt=wts, x=scaledx)/1.34)
    }
    ans <- .C(Cfunction, n = as.integer(n), m = as.integer(m), 
              r = as.integer(r), blockid=as.integer(u), 
              mu = as.double(mu[iter, , ]), 
              sigma = as.double(sigma[iter, , ]),
              x = as.double(x),
              bw = as.double(h), z = as.double(z.tmp), f = double(n*m), 
              PACKAGE="mixtools")
    lambda.f <- sweep(matrix(ans$f, ncol=m), 2, lambda[iter, ], "*")

    ## E-step (for next iteration)
    z.hat <- lambda.f/rowSums(lambda.f)
    
    loglik <- c(loglik,sum(log(rowSums(lambda.f)))) # log-likelihood
    finished <- iter >= maxiter
    if (stochastic) {
      sumpost <- sumpost + z.hat
    } 
    else if (iter > 1) { # This convergence criterion may be too simplistic:
      maxchange <- max(abs(lambda[iter,] - lambda[iter-1,]))
      if (!constbw) 
        maxchange <- max(maxchange, max(abs(bw.old - bw)))
      finished <- finished | (maxchange < eps)
    }
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, "  lambda ", round(lambda[iter, ], 4))
      cat(" time", (t1 - t0)[3], "\n")
    }
  }
  if (verb) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter, ], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }
  if (stochastic) {
    return(structure(list(data = x, posteriors = sumpost/iter, 
                          lambda = lambda, bandwidth = bw, 
                          blockid = u, lambdahat = colMeans(lambda),
                          mu = mu, muhat = colMeans(mu), 
                          sigma = sigma, sigmahat = colMeans(sigma),
                          loglik = loglik), 
                        class="spEM"))
  }
  else {
    return(structure(list(data = x, posteriors = z.hat, 
                          lambda = lambda[1:iter,], bandwidth = bw, 
                          blockid = u, lambdahat = lambda[iter,],
                          mu = mu[1:iter, , ], muhat = mu[iter, , ],
                          sigma = sigma[1:iter, , ], sigmahat = sigma[iter, , ],
                          loglik = loglik),
                        class="spEM"))
  }
}



