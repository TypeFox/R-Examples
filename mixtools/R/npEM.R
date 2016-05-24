## EM-like algorithm for a nonparametric mixture model with
## independent repeated measures - some ID, some not
npEMindrep <- # npEMindrep is an alias (only for backward compatibility)
npEM <- function(x, mu0, blockid = 1:ncol(x),
                       bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
                       h=bw, eps=1e-8,
                       maxiter=500, stochastic = FALSE, verb = TRUE){
  bw <- h # h is alternative bandwidth argument, for backward compatibility
  x <- as.matrix(x)
  n <- nrow(x)      # number of subjects
  r <- ncol(x)      # number of measurements in each subject
  u <- match(blockid, unique(blockid))
  if (is.matrix(mu0)) 
    m <- dim(mu0)[1]  # mu0=centers
  else 
    m <- mu0  # mu0=number of clusters
  if(!samebw && !is.matrix(bw)) { 
    bw <- matrix(bw, nrow=max(u), ncol=m)
  }
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
  loglik <- NULL
  orderx <- xx <- list()
  for(k in 1:max(u)) {
    xx[[k]] <- as.vector(x[, u==k])
    if (!samebw) {
      orderx[[k]] = order(xx[[k]]) # only needed for IQR calculation for bw
    }
  }
  Cfunction <- ifelse(samebw, "KDErepeated", "KDErepeatedbw")
  
  while (!finished) {
    iter <- iter + 1
    bw.old <- bw
    t0 <- proc.time()

    ## Note:  Enter loop assuming E-step is done -- i.e., z.hat in place    
    ## M-Step
    lambda[iter, ] <- colMeans(z.hat)

    ## density estimation step
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
    fkernel <- matrix(1, n, m)
    h <- bw # This is for samebw == TRUE
    for (k in 1:max(u)) {
      r2 <- sum(u == k)	# block size
      if (!samebw) {
        wts <- apply(z.tmp, 2, function(z) rep(z/r2, r2))
        variances <- colSums(wts * outer(xx[[k]], colSums(wts * xx[[k]]), '-')^2)
        iqr <- apply(as.matrix(wts[orderx[[k]],]), 2, wIQR, xx[[k]][orderx[[k]]], 
                     already.sorted=TRUE, already.normalized=TRUE)
        h <- bw[k, ] <- 0.9 * pmin(sqrt(variances), iqr/1.34) * 
                   pmax(1,r2*n*lambda[iter, ])^(-1/5) 
                   # Note:  Doesn't allow "sample size" < 1.
#      browser()
      }
      ans <- .C(Cfunction, n = as.integer(n), m = as.integer(m), 
                r = as.integer(r2), x = as.double(x[,u==k]),
                bw = as.double(h), z = as.double(z.tmp), f = double(n*m), 
              PACKAGE="mixtools")
      fkernel <- fkernel * matrix(ans$f, ncol = m)
    }
    lambda.f <- sweep(fkernel, 2, lambda[iter, ], "*")
    
    ## E-step (for next iteration)
    z.hat <- lambda.f/rowSums(lambda.f)

    loglik <- c(loglik,sum(log(rowSums(lambda.f)))) # log-likelihood
    finished <- iter >= maxiter
    if (stochastic) {
      sumpost <- sumpost + z.hat
    } 
    else if (iter > 1) { # This convergence criterion may be too simplistic:
      maxchange <- max(abs(lambda[iter,] - lambda[iter-1,]))
      if (!samebw) 
        maxchange <- max(maxchange, max(abs(bw.old - bw)))
      finished <- finished | (maxchange < eps)
    }
    if (verb) {
      t1 <- proc.time()
      cat("iteration", iter, "  lambda ", round(lambda[iter, ], 4))
      cat(" time", (t1 - t0)[3], "\n")
    }
  }
  if (!samebw) {
    rownames(bw) <- paste("block", 1:max(u))
    colnames(bw) <- paste("component", 1:m)
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
                          loglik = loglik), 
                        class="npEM"))
  }
  else {
    return(structure(list(data = x, posteriors = z.hat, 
                          lambda = lambda[1:iter,], bandwidth = bw, 
                          blockid = u, lambdahat = lambda[iter,],
                          loglik = loglik),
                           class="npEM"))
  }
}

# Not sure whether the following function is really necessary:
npEMindrepbw <- function (x, mu0, blockid = 1:ncol(x), bw =
                          bw.nrd0(as.vector(as.matrix(x))), eps = 1e-08, 
                          maxiter = 500, stochastic = FALSE, verb = TRUE){
  npEM(x=x, mu0=mu0, blockid=blockid, bw=bw, samebw=FALSE, eps=eps,
       maxiter=maxiter, stochastic=stochastic, verb=verb)
}


