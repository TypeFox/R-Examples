graf <-
function (y, x, error = NULL, weights = NULL, prior = NULL, l = NULL, opt.l = FALSE,
			theta.prior.pars = c(log(10), 1), hessian = FALSE, opt.control = list(),
			verbose = FALSE, method = c('Laplace', 'EP')) {
			
  if (opt.l) {
	# call graf recursively to optimise the lengthscale parameters
	# if l is specified, use this as the starting point
 
    nlposterior <- function(theta) {
  	  it <<- it + 1
  	  if (verbose) cat(paste('\nlengthscale optimisation iteration', it, '\n'))
  	# define calculate the negative log posterior
       # if (any(theta > theta.limit)) return(.Machine$double.xmax)
  	  l <- rep(NA, k)
  	  if (length(notfacs) > 0) l[notfacs] <- exp(theta)
  	  if (length(facs) > 0) l[facs] <- 0.01
  	  if (any(is.na(l))) stop('missing lengthscales')
  	  
        llik <- -graf(y, x, error, weights, prior = prior, l = l,
  					verbose = verbose, method = method)$mnll
        lpri <- theta.prior(theta)
        lpost <- llik + lpri
  	  if (verbose) cat(paste('\nlog posterior:', lpost, '\n'))
        lpost <- ifelse(is.finite(lpost), lpost, -.Machine$double.xmax)
        return(-lpost)
    }
	
  	k <- ncol(x)
  	# set up initial lengthscales
  	if (is.null(l)) l <- rep(1, k)
  	else if (length(l) != k) stop(paste('l must have', k, 'elements'))
  	
  	# find factors and drop them from theta
  	notfacs <- 1:k
  	facs <- which(unlist(lapply(x, is.factor)))
  	if (length(facs) > 0) {
  		notfacs <- notfacs[-facs]
  		l[facs] <- 0.01
  	}
  	theta <- log(l[notfacs])
  	
  	# if we want the hessian (for later MC integration) turn off the limit to theta 
  	#if (hessian) theta.limit = Inf
  	
  	# use hyperprior parameters
  	theta.prior <- function(theta) {
  			sum(dnorm(theta, theta.prior.pars[1], theta.prior.pars[2], log = TRUE))
  	}
  	
  	it <- 0
  	if (length(notfacs) == 1)  {
  		meth <- 'Brent'
  		low <- -100
  		up <- 100
  	} else {
  		meth <- 'BFGS'
  		low <- -Inf 
  		up <- Inf
  	}
  	# run numerical optimisation on the hyperparameters
  	# if (is.null(opt.tol)) opt.tol <- sqrt(.Machine$double.eps)
      opt <- optim(theta, nlposterior, hessian = hessian, lower = low, upper = up,
  					method = meth, control = opt.control)
  	
  	# get the resultant lengthscales
  	l[notfacs] <- exp(opt$par)
  	
  	# replace hessian with the hessian matrix or NULL
  	if(hessian) hessian <- opt$hessian
  	else hessian <- NULL
  	
  	# fit the final model and return
  	return (graf(y, x, error, weights, prior, l = l, verbose = verbose, hessian = hessian, method = method))
  } # end opt.l if statement
  
  method = match.arg(method)
  
  if (!is.data.frame(x)) stop ("x must be a dataframe")
  
  # convert any ints to numerics
  for(i in 1:ncol(x)) if (is.integer(x[, i])) x[, i] <- as.numeric(x[, i])
  
  obsx <- x
  k <- ncol(x)
  n <- length(y)
  
  if (is.null(weights)) {
    # if weights aren't provided
    weights <- rep(1, n)
  } else {
    # if they are, run some checks
    # throw an error if weights are specified with EP
    if (method == 'EP') {
      stop ('weights are not implemented for the EP algorithm (yet)')
    }
    # or if any are negative
    if (any(weights < 0)) {
      stop ('weights must be positive or zero')
    }
  }
  
  # find factors and convert them to numerics
  notfacs <- 1:k
  facs <- which(unlist(lapply(x, is.factor)))
  if (length(facs) > 0) notfacs <- notfacs[-facs]
  for (fac in facs) {
    x[, fac] <- as.numeric(x[, fac])
  }
  x <- as.matrix(x)
  
  # scale the matrix, retaining scaling
  scaling <- apply(as.matrix(x[, notfacs]), 2, function(x) c(mean(x), sd(x)))
  for (i in 1:length(notfacs)) {
    x[, notfacs[i]] <- (x[, notfacs[i]] - scaling[1, i]) / scaling[2, i]
  }

  # set up the default prior, if not specified
  exp.prev <- sum(weights[y == 1]) / sum(weights)
  if (is.null(prior))  mnfun <- function(x) rep(exp.prev, nrow(x))
  else mnfun <- prior
  
  # give an approximation to l, if not specified (or optimised)
  if (is.null(l)) {
	l <- rep(0.01, k)
	l[notfacs] <- apply(x[y == 1, notfacs, drop = FALSE], 2, sd) * 8
  }
  # calculate mean (on unscaled data and probability scale)
  mn <- mnfun(obsx)
  
  # fit model
  if (method == 'Laplace') {
    # by Laplace approximation
    fit <- graf.fit.laplace(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = verbose)
  } else {
    # or using the expectation-propagation algorithm
    fit <- graf.fit.ep(y = y, x = as.matrix(x), mn = mn, l = l, wt = weights, e = error, verbose = FALSE)
  }
  
  fit$mnfun = mnfun
  fit$obsx <- obsx
  fit$facs <- facs
  fit$hessian <- hessian
  fit$scaling <- scaling
  fit$peak = obsx[which(fit$MAP == max(fit$MAP))[1], ]
  class(fit) <- "graf"
  fit
}