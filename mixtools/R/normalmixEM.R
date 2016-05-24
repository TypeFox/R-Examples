## Use an ECM algorithm (in the sense of Meng and Rubin, Biometrika 1993)
## to search for a local maximum of the likelihood surface for a 
## univariate finite mixture of normals with possible equality
## constraints on the mean and stdev parameters.
normalmixEM <-
function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, 
          mean.constr = NULL, sd.constr = NULL,
          epsilon = 1e-08, maxit = 1000, maxrestarts=20, 
          verb = FALSE, fast=FALSE, ECM = FALSE,
          arbmean = TRUE, arbvar = TRUE) {
  warn <- options(warn=-1) # Turn off warnings
  x <- as.vector(x)
  tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
                        k = k, arbmean = arbmean, arbvar = arbvar)
  lambda <- tmp$lambda 
  mu <- tmp$mu 
  sigma <- tmp$s
  k <- tmp$k
  arbvar <- tmp$arbvar 
  arbmean <- tmp$arbmean
  if (fast==TRUE && k==2 && arbmean==TRUE) {
    a <- normalmixEM2comp (x, lambda=lambda[1], mu=mu, sigsqrd=sigma^2, 
                           eps=epsilon, maxit=maxit, verb=verb)
  } else {
    z <- parse.constraints(mean.constr, k=k, allsame=!arbmean)
    meancat <- z$category; meanalpha <- z$alpha
    z <- parse.constraints(sd.constr, k=k, allsame=!arbvar)
    sdcat <- z$category; sdalpha <- z$alpha
    ECM <- ECM || any(meancat != 1:k) || any(sdcat != 1)
    n <- length(x)
    notdone <- TRUE
    restarts <- 0
    while(notdone) {
      # Initialize everything
      notdone <- FALSE
      tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
                            k = k, arbmean = arbmean, arbvar = arbvar)
      lambda <- tmp$lambda
      mu <- tmp$mu
      k <- tmp$k
      sigma <- tmp$s
      var <- sigma^2
      diff <- epsilon+1
      iter <- 0
      postprobs <- matrix(nrow = n, ncol = k)
      mu <- rep(mu, k)[1:k]
      sigma <- rep(sigma,k)[1:k]
      # Initialization E-step here:
      z <- .C("normpost", as.integer(n), as.integer(k),
              as.double(x), as.double(mu), 
              as.double(sigma), as.double(lambda),
              res2 = double(n*k), double(3*k), post = double(n*k),
              loglik = double(1), PACKAGE = "mixtools")
      postprobs <- matrix(z$post, nrow=n)
      res <- matrix(z$res2, nrow=n)
      ll <- obsloglik <- z$loglik
      while (diff > epsilon && iter < maxit) {
        # ECM loop, 1st M-step: condition on sigma, update lambda and mu
        lambda <- colMeans(postprobs)
        mu[meancat==0] <- meanalpha[meancat==0]
        if (max(meancat)>0) {
          for(i in 1:max(meancat)) {
            w <- which(meancat==i)
            if (length(w)==1) {
              mu[w] <- sum(postprobs[,w]*x) / (n*lambda[w])
            } else {
              tmp <- t(postprobs[,w])*(meanalpha[w]/sigma[w]^2)
              mu[w] <- meanalpha[w] * sum(t(tmp)*x) / sum(tmp*meanalpha[w])
            }
          }
        }

        if (ECM) {  # If ECM==FALSE, then this is a true EM algorithm and
                    # so we omit the E-step between the mu and sigma updates
          # E-step number one:
          z <- .C("normpost", as.integer(n), as.integer(k),
                  as.double(x), as.double(mu), 
                  as.double(sigma), as.double(lambda),
                  res2 = double(n*k), double(3*k), post = double(n*k),
                  loglik = double(1), PACKAGE = "mixtools")
          postprobs <- matrix(z$post, nrow=n)
          res <- matrix(z$res2, nrow=n)

        # ECM loop, 2nd M-step: condition on mu, update lambda and sigma
          lambda <- colMeans(postprobs) # Redundant if ECM==FALSE
        }
        sigma[sdcat==0] <- sdalpha[sdcat==0]
        if (max(sdcat)>0) {
          for(i in 1:max(sdcat)) {
            w <- which(sdcat==i)
            if (length(w)==1) {
              sigma[w] <- sqrt(sum(postprobs[,w]*res[,w]) / (n*lambda[w]))
            } else {
              tmp <- t(postprobs[,w]) / sdalpha[w]
              sigma[w] <- sdalpha[w] * sqrt(sum(t(tmp) * res[,w])/ (n * sum(lambda[w])))
            }
          }
          if(any(sigma < 1e-08)) {
            notdone <- TRUE
            cat("One of the variances is going to zero; ",
                "trying new starting values.\n")
            restarts <- restarts + 1
            lambda <- mu <- sigma <- NULL
            if(restarts>maxrestarts) { stop("Too many tries!") }
            break
          }
        }
        
        # E-step number two:
        z <- .C("normpost", as.integer(n), as.integer(k),
                as.double(x), as.double(mu), 
                as.double(sigma), as.double(lambda),
                res2 = double(n*k), double(3*k), post = double(n*k),
                loglik = double(1), PACKAGE = "mixtools")
        postprobs <- matrix(z$post, nrow=n)
        res <- matrix(z$res2, nrow=n)
        newobsloglik <- z$loglik
        diff <- newobsloglik - obsloglik
        obsloglik <- newobsloglik
        ll <- c(ll, obsloglik)
        iter <- iter + 1
        if (verb) {
          cat("iteration =", iter, " log-lik diff =", diff, " log-lik =", 
              obsloglik, "\n")
          print(rbind(lambda, mu, sigma))
        }
      }
    }
    if (iter == maxit) {
      cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of iterations=", iter, "\n")
    if(arbmean == FALSE){
      scale.order = order(sigma)
      sigma.min = min(sigma)
      postprobs = postprobs[,scale.order]
      colnames(postprobs) <- c(paste("comp", ".", 1:k, sep = ""))
      a=list(x=x, lambda = lambda[scale.order], mu = mu, sigma = sigma.min, 
             scale = sigma[scale.order]/sigma.min, loglik = obsloglik, 
             posterior = postprobs, all.loglik=ll, restarts=restarts, 
             ft="normalmixEM")
    } else {
      colnames(postprobs) <- c(paste("comp", ".", 1:k, sep = ""))
      a=list(x=x, lambda = lambda, mu = mu, sigma = sigma, loglik = obsloglik, 
             posterior = postprobs, all.loglik=ll, restarts=restarts, 
             ft="normalmixEM")
    }
  }
  class(a) = "mixEM"
  options(warn) # Reset warnings to original value
  a
}

