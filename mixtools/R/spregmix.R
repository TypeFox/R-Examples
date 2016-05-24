spregmix = function (lmformula, bw = NULL, constbw = FALSE,
                     bwmult = 0.9, # gives Silverman's rule; or use 1.06 for Scott's
                     z.hat = NULL, symm = TRUE, betamethod="LS",
                     m = ifelse(is.null(z.hat), 2, ncol(z.hat)),
                     epsilon = 1e-04, maxit = 1000, verbose = FALSE, 
                     ...) 
{
# m, not k, is the number of components here
  y <- model.frame(lmformula, ...)[,1]
  x <- model.matrix(lmformula, ...)
  n <- length(y)
  p <- ncol(x)
  tt0 <-  proc.time() # for total time

  # Now initialize z.hat, which is the only thing needed to start the iterations
  # as long as we use a method such as least squares to find beta
  if(is.null(z.hat)) {
    z.hat <- matrix(runif(n*m), n, m)
    z.hat <- sweep(z.hat, 1, rowSums(z.hat), "/")
  }

  L1norm <- function(beta, y, x, p) sum(p*abs(y - x %*% beta))
  nploglik <- function(beta, y, x, p, bw, res, wts, symm) {
    sum(p*log(wkde(x=res, u=as.vector(y-x%*%beta),  w=wts, bw=bw, sym=symm))) 
  }
  lambda <- matrix(0, maxit, m)
  beta <- matrix(0, p, m)
  for (j in 1:m) { # Initialize beta to weighted LS solution
        wx <- sweep(x, 1, z.hat[,j], "*")
        beta[,j] <- solve(t(wx) %*% x, t(wx) %*% y) 
  }
  iter <- 0
  finished <- FALSE
#  j <- 0
#plot(tonedata)
  SPweight <- ifelse(betamethod=="transition", 0, 1)
  while (!finished) {
#    j <- j + 1; if (j>m) {
#browser();
#plot(x[,2],y);
#      j <- 1}
    iter <- iter + 1
    t0 <- proc.time()
    ## M-step
    lambda[iter, ] <- colMeans(z.hat)
    oldbeta <- beta
    for (j in 1:m) {
#abline(beta[1,j],beta[2,j], col=1+j, lwd=3, lty=2)
      if (betamethod=="LS") {
        wx <- sweep(x, 1, z.hat[,j], "*")
        beta[,j] <- solve(t(wx) %*% x, t(wx) %*% y) # Weighted least squares solution
      } else if (betamethod=="L1") { # Weighted L1 solution
        beta[,j] <- optim(par=beta[,j], fn=L1norm, y=y, x=x, p=z.hat[,j])$par
      } else if (betamethod=="transition") { # transition from LS to NP
        wx <- sweep(x, 1, z.hat[,j], "*")
        beta[,j] <- SPweight * 
                    optim(par=beta[,j], fn=nploglik, y=y, x=x, p=z.hat[,j], 
                          bw=bw, res=as.vector(y-x%*%beta), 
                          wts=as.vector(z.hat), symm=symm, 
                          control=list(fnscale=-1))$par + # NP loglik
                    (1-SPweight) *
                    solve(t(wx) %*% x, t(wx) %*% y) # Weighted least squares 
        SPweight <- min(1, SPweight + (1e-4)*2^(iter-1))
      } else { # Nonparametric loglikelihood 
        beta[,j] <- optim(par=beta[,j], fn=nploglik, y=y, x=x, p=z.hat[,j], 
                          bw=bw, res=as.vector(y-x%*%beta), 
                          wts=as.vector(z.hat), symm=symm, 
                          control=list(fnscale=-1))$par
#res=as.vector(y-x%*%beta)
#plot(res,wkde(x=res, u=res,  w=as.vector(z.hat), bw=bw, sym=symm))
#browser()

      }
      ## Here, we might try other methods of estimating beta.
    }

    ## update the bandwidth, if necssary
    xbetavec <- as.double(x %*% beta)
    zhatvec <- as.double(z.hat)
    if (is.null(bw) || !constbw) {
      res <- y-xbetavec
      np.sigma <- sqrt(sum(res^2 * zhatvec)/(n-1))
      bw <- bwmult / n^(1/5) * min(np.sigma, wiqr<-wIQR(wt=zhatvec, x=res)/1.34) 
#print(c(np.sigma, wiqr))
    }

    ## density estimation step  
    if (symm) {
      ans <- .C ("KDEsymloc2", n=as.integer(n), m=as.integer(m),
                 mu = xbetavec, y=as.double(y), bw=as.double(bw),
                 z=as.double(z.hat), f=double(n*m))
    } else {
      ans <- .C ("KDEloc2", n=as.integer(n), m=as.integer(m),
                 mu = xbetavec, y=as.double(y), bw=as.double(bw),
                 z=as.double(z.hat), f=double(n*m))
    }
    fkernel <- matrix(ans$f, ncol=m)
    lambda.f <- sweep(fkernel, 2, lambda[iter,], "*")

    ## E-step (for next iteration)
    z.hat <- lambda.f/rowSums(lambda.f)    
    
    ## Test convergence criteria
    finished <- (iter > 1 && SPweight == 1 &&
                 (iter >= maxit  ||
                  max(max(abs(beta-oldbeta)),
                      max(abs(lambda[iter,]-lambda[iter-1,]))) < epsilon
                  ))
    
    ## Print message for each iteration
    if (verbose) {
      t1 <- proc.time()
      cat("iteration ", iter, "  lambda ", round(lambda[iter,], 4), 
          "  bandwidth ", round(bw, 3))
      cat(" time", (t1 - t0)[3], "\n")
    }
  }

  ## Print final message after convergence
  if (verbose) {
    tt1 <- proc.time()
    cat("lambda ", round(lambda[iter,], 4))
    cat(", total time", (tt1 - tt0)[3], "s\n")
  }

  res <- ans$y - ans$mu
  xx <- res
  yy <- ans$f
  if (symm) {
    xx <- c(xx, -res) # exploiting symmetry
    yy <- c(yy, yy)
  }
  ox <- order(xx)
  np.sigma <- sqrt(sum(res^2 * ans$z)/(n-1))
  
  loglik <- rep(0, n)
  for(j in 1:m) {
    loglik <- loglik + z.hat[,j]*wkde(as.vector(y-x%*%beta[,j]),
                                       w=as.vector(z.hat[,j]), bw=bw, sym=symm)
  }
  loglik <- sum(log(loglik))
  
  a <- list(x=x, y=y, lambda = lambda[1:iter,], beta = beta,
         posterior = z.hat, np.stdev = np.sigma, bandwidth=bw,
         density.x = xx[ox], density.y = yy[ox], symmetric=symm,
         loglik = loglik,
         ft="regmixEM")
  class(a) = "npEM"
  a
} 
