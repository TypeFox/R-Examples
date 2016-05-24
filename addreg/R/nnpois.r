nnpois <- function(y, x, standard, offset, start, control = list()) 
{
  control <- do.call("addreg.control", control)
  x <- as.matrix(x)
  xnames <- dimnames(x)[[2L]]
  ynames <- if(is.matrix(y))
    rownames(y)
  else names(y)
  
  if(any(x < 0)) stop("x must be non-negative")
  if(any(apply(x, 2, function(col) all(col==0)))) stop("x contains column with all 0")
  
  nvars <- ncol(x)
  nobs <- NROW(y)
  
  fam <- poisson(link = identity)
  eval(fam$initialize)
  
  mu.eta <- fam$mu.eta
  linkinv <- fam$linkinv
  dev.resids <- fam$dev.resids
  aic <- fam$aic
  loglik <- function(y, mu) sum(y*log(mu) - mu - lgamma(y+1))
  
  weights <- rep(1, nobs)
  if (is.null(standard)) standard <- rep.int(1, nobs)
  if (is.null(offset)) offset <- rep.int(0, nobs)
  if (any(offset < 0))
    stop("offset must be non-negative")
  
  converged <- FALSE
  
  coefold <- if(!is.null(start)) {
    if (length(start) != nvars)
        stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                nvars, paste(deparse(xnames), collapse = ", ")),
            domain = NA)
    else if(any(start <= control$bound.tol))
        stop("'start' is on our outside the boundary of the parameter space (consider 'bound.tol')", domain = NA)
	else start
    } else {
        simple <- mean(y / standard) / colMeans(x) + 2*control$bound.tol
        trymat <- tryCatch(as.vector(solve(t(x)%*%x) %*% t(x) %*% (y)) + 2*control$bound.tol, error = function(e) NULL)
        if(is.null(trymat)) simple
        else if(any(trymat < control$bound.tol)) simple
        else trymat
    }

  eta <- drop(x %*% coefold) + offset
  mu <- standard * linkinv(eta)
  dev.old <- sum(dev.resids(y, mu, weights))
  
  if(control$trace) cat("Deviance =", dev.old, "Iterations - 0\n")
  
  std.div <- 1 / colSums(standard * x)
  
  for(iter in 1L:control$maxit) {
    
    y.over.fits <- y / linkinv(eta)
    y.over.fits[linkinv(eta)==0] <- 0
    
    coefnew <- coefold * colSums(y.over.fits * x) * std.div
    
    eta <- drop(x %*% coefnew) + offset
    mu <- standard * linkinv(eta)
    dev.new <- sum(dev.resids(y, mu, weights))
    
    ll.new <- loglik(y, mu)
    
    if(control$trace) cat("Deviance =", dev.new, "Iterations -", iter, "\n")
    
    if(conv.test(coefold, coefnew, control$epsilon)) {
      converged = TRUE
      break
    }
    
    coefold <- coefnew
    dev.old <- dev.new
  }
  
  residuals <- (y-mu)/mu.eta(eta)
  
  names(coefnew) <- xnames
  names(residuals) <- ynames
  names(mu) <- ynames
  names(eta) <- ynames
  names(y) <- ynames
  
  aic.model <- aic(y, nobs, mu, weights, dev.new) + 2 * nvars
  aic.c.model <- aic.model + 2 * nvars * (nvars + 1) / (nobs - nvars - 1)
  
  wtdmu <- standard * sum(weights * y / standard)/sum(weights)
  nulldev <- sum(dev.resids(y, wtdmu, weights))
  nulldf <- nobs - 1
  resdf <- nobs - nvars
  
  boundary <- any(coefnew < control$bound.tol)
  
  list(coefficients = coefnew, residuals = residuals, fitted.values = mu, rank = nvars, family = fam,
       linear.predictors = eta, deviance = dev.new, aic = aic.model, aic.c = aic.c.model,
       null.deviance = nulldev, iter = iter, weights = weights, prior.weights = weights,
       standard = standard, df.residual = resdf, df.null = nulldf, y = y, 
       converged = converged, boundary = boundary, loglik = ll.new, nn.design = x)
  
}