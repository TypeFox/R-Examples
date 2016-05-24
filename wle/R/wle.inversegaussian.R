#############################################################
#                                                           #
#	wle.inversegaussian function                        #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 18, 2011                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

wle.inversegaussian <- function(x, boot=30, group, num.sol=1, raf="HD", smooth=0.003, tol=10^(-6), equal=10^(-3), max.iter=500, use.smooth=TRUE, tol.int, verbose=FALSE) {

  x <- as.vector(x)
  size <- length(x)
  result <- list()  
  mufun <- function(x, w=rep(1, length(x))) {
    result <- c(w%*%x)/sum(w)
    return(result)
  }
  lambdafun <- function(x, mu, w=rep(1, length(x))) {
    result <- mu^2*sum(w)/sum(w*((x-mu)^2)/x)
    return(result)
  }
  
  raf <- switch(raf,
	HD = 1,
	NED = 2,
	SCHI2 = 3,
	-1)

  if (raf==-1)
    stop("Please, choose the RAF: HD=Hellinger Disparity, NED=Negative Exponential Disparity, SCHI2=Symmetric Chi-squares Disparity")

  if (missing(group))
    group <- 0

  if (size<2)
    stop("Number of observation must be at least equal to 2")

  if (group<2) {
    group <- max(round(size/8),2)
    if (verbose) cat("wle.inversegaussian: dimension of the subsample set to default value: ",group,"\n")
  }

  maxboot <- sum(log(1:size))-(sum(log(1:group))+sum(log(1:(size-group))))

  if (boot<1 | log(boot) > maxboot)
    stop("Bootstrap replication not in the range")

  if (!(num.sol>=1)) {
    if (verbose) cat("wle.inversegaussian: number of solution to report set to 1 \n")
    num.sol <- 1
  }

  if (max.iter<1) {
    if (verbose) cat("wle.inversegaussian: max number of iteration set to 500 \n")
    max.iter <- 500
  }

  if (smooth<10^(-5))
    if (verbose) cat("wle.inversegaussian: the smooth parameter seems too small \n")

  if (tol<=0) {
    if (verbose) cat("wle.inversegaussian: the accuracy must be positive, using default value: 10^(-6) \n")
    tol <- 10^(-6)
  }

  if (equal<0) {
    if (verbose) cat("wle.inversegaussian: the equal parameter must be greater than tol, using default value: tol+10^(-3) \n")
    equal <- tol+10^(-3)
  }

  if (!is.logical(use.smooth)) {
    if (verbose) cat("wle.inversegaussian: the use.smooth must be a logical value, using default value \n")
    use.smooth <- TRUE
  }

  if (missing(tol.int)) {
    tol.int <- tol*10^(-4)
  } else {
    if (tol.int <=0) {
      if (verbose) cat("wle.inversegaussian: tol.int must be positive, using default value \n")
       tol.int <- tol*10^(-4) 
    } 
  }
  
  tot.sol <- 0
  not.conv <- 0
  iboot <- 0

  while (tot.sol < num.sol & iboot < boot) {
    iboot <- iboot + 1
    xboot <- x[sample(size, group)]
    mu <- mufun(xboot)
    lambda <- lambdafun(xboot, mu)

    if (lambda > 10000 & use.smooth) {
      lambda <- 1000
      if (verbose)
        cat("wle.inversegaussian: some initial values for 'lambda' are greater than 10000 reset them to 1000. Try use.smooth=FALSE if you do not want them to be truncated. \n")
    }
      
    
    xdiff <- tol + 1
    iter <- 0
    while (xdiff > tol & iter < max.iter) {
      iter <- iter + 1
#      cat('mu ', mu, '\n')
#      cat('lambda ', lambda, '\n')      
      muold <- mu
      lambdaold <- lambda      
      dsup <- max(x)+ 3*smooth/lambda      
      z <- .Fortran("wleinvga",
	    as.double(x),
            as.double(x),
	    as.integer(size),
	    as.integer(size),                     
	    as.integer(raf),
            as.double(1),
	    as.double(smooth/lambda),
            as.integer(1*use.smooth),
            as.double(dsup),
	    as.double(tol),
            as.double(tol.int),
	    as.double(mu),
	    as.double(lambda),
	    weights=double(size),
	    density=double(size),
	    model=double(size),
            PACKAGE = "wle")
      mu <- mufun(x, z$weights)
      lambda <- lambdafun(x, mu, z$weights)      
      xdiff <- max(abs(mu - muold), abs(lambda - lambdaold))
    }
#### end of while (xdiff > tol & iter < max.iter)

    if (iter < max.iter) {
      if (tot.sol==0) {
        mu.store <- mu
        l.store <- lambda
        w.store <- z$weights
        t.store <- sum(z$weights)/size
        f.store <- z$density
        m.store <- z$model
        d.store <- f.store/m.store-1
        tot.sol <- 1
      } else {
        if (min(abs(m.store-mu), abs(l.store-lambda)) > equal) {
          mu.store <- c(mu.store, mu)
          l.store <- c(l.store, lambda)
          w.store <- rbind(w.store, z$weights)
          t.store <- c(t.store, sum(z$weights)/size)
          f.store <- rbind(f.store, z$density)
          m.store <- rbind(m.store, z$model)
          d.store <- rbind(d.store, z$density/z$model - 1)
          tot.sol <- tot.sol + 1
        }
      }

    } else not.conv <- not.conv + 1  
  }
##### end of while (tot.sol < num.sol & iboot < boot)

  if (tot.sol) {
    result$mu <- c(mu.store)
    result$lambda <- c(l.store)
    result$tot.weights <- c(t.store)
    result$weights <- w.store
    result$f.density <- f.store
    result$m.density <- m.store
    result$delta <- d.store
    result$tot.sol <- tot.sol
    result$not.conv <- not.conv
  } else {
    if (verbose)
      cat("wle.poisson: No solutions are fuond, checks the parameters\n")
    result$mu <- NA
    result$lambda <- NA
    result$tot.weights <- NA
    result$weights <- rep(NA,size)
    result$f.density <- rep(NA,size)
    result$m.density <- rep(NA,size)
    result$delta <- rep(NA,size)
    result$tot.sol <- 0
    result$not.conv <- boot
  }

  result$call <- match.call()
  class(result) <- "wle.inversegaussian"
  return(result)
}

#############################################################
#                                                           #
#	print.wle.inversegaussian function                  #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 18, 2011                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2011 Claudio Agostinelli              #
#                                                           #
#############################################################

print.wle.inversegaussian <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Location:\n")
    print.default(format(x$mu, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Precision:\n")
    print.default(format(x$lambda, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nNumber of solutions ", x$tot.sol, "\n")
    cat("\n")
    invisible(x)
}



