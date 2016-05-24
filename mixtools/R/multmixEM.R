multmixEM <- function (y, lambda = NULL, theta = NULL, k = 2, maxit = 10000, 
    epsilon = 1e-08, verb = FALSE) {
  if (class(y)=="list" && !is.null(y$y)) {
    y <- y$y
  }
  n <- nrow(y)
  p <- ncol(y)
  m <- colSums(y)
  r <- rowSums(y) # These need not all be the same

  tmp <- multmix.init(y=y, lambda=lambda, theta=theta, k=k)
  lambda <- tmp$lambda
  theta <- tmp$theta
  k <- tmp$k
  restarts<-0
  mustrestart <- FALSE

  llconstant <- sum(lgamma(r+1)) - sum(lgamma(y+1))
  while (restarts < 50) {
    ll <- NULL
    iter <- 0
    diff <- epsilon+1 # to ensure entering main EM loop
    # Normalize rows of theta matrix
    theta <- theta/rowSums(theta)
    theta <- pmax(theta, 1e-100) # ensure no zeros
    # preparatory E-step prior to entering main EM loop
    loglamcd <- log(lambda) + log(theta) %*% t(y) # kxn matrix of log(lambda * component densities)
    z <- .C("multinompost", as.integer(n), as.integer(k),
        as.double(loglamcd), post=double(n*k), loglik=as.double(llconstant),
        PACKAGE = "mixtools")
    post <- matrix(z$post, ncol=k)
    newll <- z$loglik
	tmp.post <- (post==0)
	if(any(apply(tmp.post,2,sum)==n)){
		diff <- epsilon
		mustrestart <- TRUE
	}
    while ((iter < maxit) && diff > epsilon) {  # main EM loop
      iter <- iter + 1
      oldll <- newll
      ll <- c(ll, oldll)
      # M-step:  First update theta values (proportional to sum of posteriors * data)
      theta <- t(post) %*% y
      theta <- theta/rowSums(theta)
      theta <- pmax(theta, 1e-100) # ensure no zeros
      # M-step:  Update the lambdas as usual for a finite mixture
      lambda <- colMeans(post)      
      # E-step:  prepare to find posteriors using C function
      loglamcd <- log(lambda) + log(theta) %*% t(y) # kxn matrix of log(lambda * component densities)
      # E-step:  Call C function to return nxk matrix of posteriors along with loglikelihood
      z <- .C("multinompost", as.integer(n), as.integer(k),
              as.double(loglamcd), post=double(n*k), loglik=as.double(llconstant),
              PACKAGE = "mixtools")
      post <- matrix(z$post, ncol=k)
      newll <- z$loglik
      diff <- newll - oldll      
      if (diff<0 || is.na(newll) || is.infinite(newll) || is.nan(newll)) {
        mustrestart <- TRUE
        break
      }
      if (verb) {
        cat("iteration=", iter, "diff=", diff, "log-likelihood", 
            ll[iter], "lambda", lambda, "\n") 
      }
    }
    if (mustrestart) {
      cat("Restarting due to numerical problem.\n")
      mustrestart <- FALSE
      restarts <- restarts + 1
      tmp <- multmix.init(y=y, k=k)
      lambda <- tmp$lambda
      theta <- tmp$theta
      k <- tmp$k
    } else {
      if (iter == maxit) {
        cat("Warning: Not convergent after", maxit, "iterations\n")
      }
      theta[,p] <- 1-apply(as.matrix(theta[,1:(p-1)]),1,sum)
      colnames(theta) <- c(paste("theta", ".", 1:p, sep = ""))
      rownames(theta) <- c(paste("comp", ".", 1:k, sep = ""))
      colnames(post) <- c(paste("comp", ".", 1:k, sep = ""))
      cat("number of iterations=", iter, "\n")
      out <-list(y=y, lambda = lambda, 
                  theta = theta, loglik = ll[length(ll)], posterior = post, 
                  all.loglik=ll, restarts=restarts, ft="multmixEM")
      class(out) <- "mixEM"
      return(out)
    }
  }
  stop("Too many restarts!")
}

