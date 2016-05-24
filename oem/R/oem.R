oemfit <- function(formula, data = list(), lambda = NULL, nlambda = 100,
                   lambda.min.ratio = NULL, tolerance = 1e-3,
                   maxIter = 1000, standardized = TRUE, numGroup = 1,
                   penalty = c("lasso", "scad", "ols", "elastic.net",
                     "ngarrote", "mcp"), alpha = 3,
                   evaluate = 0, condition = -1) {
  # prepare the generic arguments
  this.call <- match.call()
  penalty <- match.arg(penalty)
  mf <- model.frame(formula = formula, data = data)
  x <- model.matrix(attr(mf, "terms"), data = mf)
  y <- model.response(mf)
  nobs <- nrow(x)
  nlambda <- as.integer(nlambda)
  tolerance <- as.double(tolerance)
  maxIter <- as.integer(maxIter)
  numGroup <- as.integer(numGroup)

  if (!standardized) {
    meanx <- apply(x, 2, mean)
    normx <- sqrt(apply((t(x) - meanx)^2, 1, sum) / nobs)
    nz <- which(normx > .0001)
    xx <- scale(x[,nz], meanx[nz], normx[nz])
    yy <- y - mean(y)
  } else {
    xx <- x
    yy <- y
  }
  nvars <- as.integer(ncol(xx))
  
  # NOTES: reset lambda.max here
  lambda.max <- max(abs(t(xx) %*% yy / nobs) * 1.1)
  if (is.null(lambda)) {
    if (is.null(lambda.min.ratio)) 
      lambda.min.ratio = ifelse(nobs < nvars, .05, 1e-3)
    if (lambda.min.ratio >= 1)
      stop("lambda.min.ratio should be less than 1")
    wlambda <- exp( seq(log(as.double(lambda.max)),
                        log(as.double(lambda.min.ratio)),
                        log(as.double(lambda.min.ratio / lambda.max))
                        /nlambda) )
    wlambda <- wlambda[1:nlambda]
  } else {
    if (any(lambda < 0)) stop("lambda should be non-negative")
    wlambda <- as.double(rev(sort(lambda)))
    nlambda <- length(wlambda)
  }
  method <- as.integer(switch(penalty,
                              ols = 0,
                              lasso = 1,
                              scad = 2,
                              elastic.net = 3,
                              ngarrote = 4,
                              mcp = 5))
  
  # determine which condition to calculate
  if (condition < 0)
    condition = ifelse(2 * nobs <= nvars, 0, 1)

  if (penalty == "ols") wlambda = 0
  result <- .Call("oemfit", xx, yy, maxIter, tolerance, wlambda,
                  method, numGroup, as.double(alpha), as.integer(evaluate),
                  as.integer(condition),
                  PACKAGE = "oem")
  result$lambda <- wlambda
  result$call <- this.call
  result$sumSquare <- apply( (yy - xx %*% result$beta)^2, 2, sum) / nobs
  
  # unstandadize
  if (!standardized && nrow(beta) == nz + 1) {
    beta <- matrix(0, ncol(xx) + 1, length(wlambda))
    beta[nz+1,] <- result$beta / normx[nz]
    beta[1,] <- mean(y) - crossprod(meanx, beta[-1,, drop = FALSE])
    result$beta <- beta
  }
  
  class(result) <- c("oemfit", class(result))
  result
}

