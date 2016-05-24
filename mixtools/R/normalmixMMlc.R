###################################################################
## "EC-MM" algorithm 
## to search for a local maximum of the likelihood surface for a 
## univariate finite mixture of normals with possible 
## linear constraints on the mean and variance parameters.
##
## (EC-MM is ECM in the sense of Meng and Rubin, 
## Biometrika 1993, where the M step is replaced by a 
## Conditional MM step)
##
## version allowing for linear constraints on the mean
## mu = M beta + C, where M is matrix(k,p) and C is vector()
## C-MM step required for the linear constraint on the Variances
## var.lincstr = matrix A (k,q) s.t.  iv = A g, where "i"nverse "v"ar
## iv = k-vector of 1/v_j's, and g = q-vector of unknown parameters
## no fast option, & ECM-MM algorithm forced (no ECM option available)
## init values for gparam are *required* here for the MM algorithm
## # default value for A could be diag(k) so that iv=g

normalmixMMlc <- function (x, 
			lambda = NULL, mu = NULL, sigma = NULL, k = 2,
			mean.constr = NULL,
	  		mean.lincstr = NULL, mean.constant = NULL,
	  		var.lincstr = NULL, gparam = NULL,
          	epsilon = 1e-08, maxit = 1000, maxrestarts=20, 
          	verb = FALSE) {
  ECM <- TRUE  # always required in this case
  A <- var.lincstr
  x <- as.vector(x);     n <- length(x)
  tmp <- normalmix.init(x = x, lambda = lambda, mu = mu, s = sigma, 
                        k = k) # no arbmean & arbvar parameters
  lambda <- tmp$lambda; mu <- tmp$mu; sigma <- tmp$s; k <- tmp$k
  # arbvar <- tmp$arbvar; arbmean <- tmp$arbmean
	arbmean <- arbvar <- TRUE # forced for parse.constraints()
	warn <- options(warn=-1) # Turn off warnings (for parsing only)
    z <- parse.constraints(mean.constr, k=k, allsame=!arbmean)    
    options(warn) # Reset warnings to original value
    meancat <- z$category; meanalpha <- z$alpha 
    
    if (!is.null(mean.lincstr)) {   # linear constraint on the means
		M <- mean.lincstr
		cat("linear constraint mu = M beta + C version\n")
		p <- dim(M)[2] # nb of columns = size of constr. mean parameter
		if (dim(M)[1] != k) 
			stop("mean.lincstr and mu dimensions must agree")
		if (is.null(mean.constant)) 
			C <- matrix(0,k,1) else C <- matrix(mean.constant,k,1)
   		}
    notdone <- TRUE
    while(notdone) {
      # Initialize everything
      notdone <- FALSE
      tmp <- normalmix.init(x=x, lambda=lambda, mu=mu, s=sigma, k=k)
      lambda <- tmp$lambda; mu <- tmp$mu; k <- tmp$k; sigma <- tmp$s

      q <- dim(A)[2] # nb of inverse variance parameters (gamma)
      g <- gparam # init g values *required* for MM algo
	  iv <- A %*% g # inverse variances, as a one-column matrix
	  v <- 1/iv
	  # is it necessary to redefined sigma from g here ?
	  sigma <- as.vector(sqrt(v))      
      diff <- epsilon+1
      iter <- 0
      postprobs <- matrix(nrow = n, ncol = k)
      restarts <- 0
      mu <- rep(mu, k)[1:k]  # is this useful? 
      sigma <- rep(sigma,k)[1:k] # sigma still needed for post computation   
      
      ## Initialization E-step here:
      z <- .C("normpost", as.integer(n), as.integer(k),
              as.double(x), as.double(mu), 
              as.double(sigma), as.double(lambda),
              res2 = double(n*k), double(3*k), post = double(n*k),
              loglik = double(1), PACKAGE = "mixtools")
      postprobs <- matrix(z$post, nrow=n)
      res <- matrix(z$res2, nrow=n) # n,m matrix of squared residuals (x_i-mu_j)^2
      ll <- obsloglik <- z$loglik
      
      ## EC-MM iterations
      while (diff > epsilon && iter < maxit) {      
        # ECM loop, 1st M-step: condition on sigma, update lambda and mu :
        lambda <- colMeans(postprobs) # update for lambda
        # update for mu, depending on constraint type
        mu[meancat==0] <- meanalpha[meancat==0]
        if (max(meancat)>0 && is.null(mean.lincstr)) { # simple constraint
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
	
		if (!is.null(mean.lincstr)) {  # linear constraint mu = M beta + C
			iv2 <- as.vector(iv)
			# A1_j = sum_i p_ij x_i/v_j
			A1 <- apply(postprobs*matrix(x,n,k),2,sum)*iv2
			B <- diag(apply(postprobs,2,sum)*iv2)
			Stemp <- solve(t(M) %*% B %*% M)
			beta <- Stemp %*% t(M) %*% (A1 - B %*% C)
			mu <- as.vector(M %*% beta + C) # coerce to vector 
			}

          # ECM E-step number one:
          z <- .C("normpost", as.integer(n), as.integer(k),
                  as.double(x), as.double(mu), 
                  as.double(sigma), as.double(lambda),
                  res2 = double(n*k), double(3*k), post = double(n*k),
                  loglik = double(1), PACKAGE = "mixtools")
          postprobs <- matrix(z$post, nrow=n)
          res <- matrix(z$res2, nrow=n)

        #### ECM loop 2nd M-step: condition on mu, update lambda 
        #### and sigma via the MM step           
        lambda <- colMeans(postprobs)          
        # Update variances with the MM algo on g, conditional on mu
        # note: code in too much steps/details for debugging
		# computing q-vector of denominators
		r0 <- postprobs*res # (i,j)th is p_{ij}*(x_i - mu_j)^2
		r1 <- r0 %*% A
		# r1(n,q) matrix, (i,l) term =\sum_j p_{ij}A_{jl}(x_i - mu_j)^2
		den <- colSums(r1) # q-vector of denominators for updating g 
		# computing q-vector of numerators
		r3 <- matrix(v,nrow=k,ncol=q)*A # (k,q) matrix of v_j.A_{jl}
		r4 <- postprobs %*% r3 # (n,q) matrix of \sum_j p_{ij} A_{jl} v_j
		num <- colSums(r4) #
		# update of gamma parameters which gives iv, v and sigma 
		g.new <- g*(num/den) 	
		iv <- A %*% g.new
		v <- 1/iv 	# needed for next computation of r3
		sigma <- as.vector(sqrt(v)) # needed in next E-steps
		g <- g.new # iterates
		                
        ## test for variance degeneration
        if(any(sigma < 1e-08)) {
            notdone <- TRUE
            cat("One of the variances is going to zero; ",
                "trying new starting values.\n")
            restarts <- restarts + 1
            lambda <- mu <- sigma <- NULL  # WHAT TO DO WITH g in this case?
            if(restarts>maxrestarts) { stop("Too many tries!") }
            break
          }        
        
        # ECM E-step number two:
        z <- .C("normpost", as.integer(n), as.integer(k),
                as.double(x), as.double(mu), 
                as.double(sigma), as.double(lambda),
                res2 = double(n*k), double(3*k), post = double(n*k),
                loglik = double(1), PACKAGE = "mixtools")
        postprobs <- matrix(z$post, nrow=n)
        res <- matrix(z$res2, nrow=n)
        newobsloglik <- z$loglik
        diff <- newobsloglik - obsloglik # does not increase that one?
        obsloglik <- newobsloglik
        ll <- c(ll, obsloglik)
        iter <- iter + 1
        if (verb) {
          cat("iteration", iter, ": log-lik diff =", round(diff,4), 
          		" log-lik =", obsloglik, "\n")
        }
      }
    }

# summurizing and returning structure
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
             gamma = g,
             ft="normalmixMMlc")
    } else {
      colnames(postprobs) <- c(paste("comp", ".", 1:k, sep = ""))
      a=list(x=x, lambda = lambda, mu = mu, sigma = sigma, loglik = obsloglik, 
             posterior = postprobs, all.loglik=ll, restarts=restarts, 
             beta = beta, gamma = g,
             ft="normalmixMMlc")
    }
  
  class(a) = "mixEM"
  a
}

