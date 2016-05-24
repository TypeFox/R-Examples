## function to set some control parameters for boosting (replaces crch.control() for boosting)
crch.boost <- function(maxit = 100, nu = 0.1, start = NULL, 
  dot = "separate", mstop = c("max", "aic", "bic", "cv"),  nfolds = 10, 
  foldid = NULL)
{
  if(is.numeric(mstop)) {
    maxit <- mstop
    mstop <- "max"
  }
  rval <- list(maxit = maxit, nu = nu, nfolds = nfolds, foldid = foldid, 
    start = start, dot = dot, mstop = match.arg(mstop), fit = "crch.boost.fit")
  rval
}


## Function to (re)standardize vectors and model matrices:
standardize.matrix <- function(x, restandardize = FALSE, center = NULL, scale = NULL,
  weights = rep(1, NROW(x))) {
  if(is.null(center)) center <- attr(x, "standardize:center")
  if(is.null(scale)) scale <- attr(x, "standardize:scale")
  x <- as.matrix(x)
  if(is.null(center))
    center <- colSums(x*weights)/sum(weights)

  center[grep("(Intercept)", colnames(x))] <- 0
  mcenter <- matrix(rep(center, nrow(x)), nrow(x), ncol(x), byrow = TRUE)

  if(is.null(scale))
    scale <- sqrt(colSums(weights*(x-mcenter)^2)/(sum(weights)))

  scale[grep("(Intercept)", colnames(x))] <- 1
  mscale  <- matrix(rep(scale , nrow(x)), nrow(x), ncol(x), byrow = TRUE)
 
  if(!restandardize) { 
    x <- (x - mcenter)/mscale
  } else {
    x <- x*mscale + mcenter
  }
  attr(x, "standardize:center") <- center
  attr(x, "standardize:scale") <- scale
  x  
}

## Function to (re)standardize coefficients
standardize.coefficients <- function(coef, restandardize = FALSE, center, scale) {
  interceptind <- grep("(Intercept)", names(coef))
  if(length(interceptind > 0)) {
    intercept <- coef[interceptind]
    center <- center[-interceptind]
    scale <- scale[-interceptind]
    coef <- coef[-interceptind]
  } else {
    intercept <- 0
    names(intercept) <- "(Intercept)"
  }
  if(restandardize) {
    intercept <- intercept - sum(coef*center/scale)
    coef <- coef/scale
  } else {
    coef <- coef*scale
    intercept <- intercept + sum(coef*center/scale)
  }
  c(intercept, coef)  
}


## actual fitting function
crch.boost.fit <- function(x, z, y, left, right, truncated = FALSE, 
  dist = "gaussian", df = NULL, link.scale = "log",
  weights = NULL, offset = NULL, control = crch.boost()) 
{
  ## response and regressor matrix
  n <- NROW(x)  
  k <- NCOL(x)
  if(is.null(weights)) weights <- rep.int(1, n)
  nobs <- sum(weights > 0)
  dfest <- identical(dist, "student") & is.null(df)
  if(is.null(offset)) offset <- rep.int(0, n)
  if(!is.list(offset)) offset <- list(location = offset, scale = rep.int(0, n))
  if(is.null(z)) {
    q <- 1L
    z <- matrix(1, ncol = q, nrow = n)
    colnames(z) <- "(Intercept)"
    rownames(z) <- rownames(x)
  } else {
    q <- NCOL(z)
    if(q < 1L) stop("scale regression needs to have at least one parameter")
  }

  ## control parameters
  maxit <- control$maxit
  nu <- control$nu
  start <- control$start
  mstop <- control$mstop

  ## extend left and right to vectors of length n
  if(length(left) == 1) left <- rep(left, n)
  if(length(right) == 1) right <- rep(right, n)

  ## link
  if(is.character(link.scale)) {
    linkstr <- link.scale
    if(linkstr != "quadratic") {
      linkobj <- make.link(linkstr)
    } else {
      linkobj <- structure(list(
        linkfun = function(mu) mu^2,
        linkinv = function(eta) sqrt(eta),
        mu.eta = function(eta) 1/2/sqrt(eta),
        valideta = function(eta) TRUE,
        name = "quadratic"
      ), class = "link-glm")
    }
  } else {
    linkobj <- link.scale
    linkstr <- link.scale$name
  }
  linkfun <- linkobj$linkfun
  linkinv <- linkobj$linkinv
  mu.eta <- linkobj$mu.eta


  ## scale response and regressors
  x <- standardize.matrix(x, weights = weights)
  z <- standardize.matrix(z, weights = weights)
  basefit <- crch.fit(x[,1, drop = FALSE], z[,1, drop = FALSE], y, left, right, 
      truncated, dist, df, link.scale, weights, offset)
  y <- standardize.matrix(y, center = basefit$coefficients$location, scale = linkinv(basefit$coefficients$scale))
    
  standardize <- list(
    x = list(center = attr(x, "standardize:center"), scale = attr(x, "standardize:scale")),
    z = list(center = attr(z, "standardize:center"), scale = attr(z, "standardize:scale")),
    y = list(center = attr(y, "standardize:center"), scale = attr(y, "standardize:scale")))

  ## standardize left, right, and offset
  left <- (left - standardize$y$center)/standardize$y$scale
  right <- (right - standardize$y$center)/standardize$y$scale
  offset[[1L]]  <- offset[[1L]]/standardize$y$scale

    
  ## restandardize coefpath if start is supplied
  if(!is.null(start)) {
    beta <- start$coefpath[,seq.int(length.out = k)]
    beta <- t(apply(beta, 1, standardize.coefficients, center = standardize$x$center, scale = standardize$x$scale))
    beta[, grep("(Intercept)", colnames(beta))] <- beta[, grep("(Intercept)", colnames(beta))] - standardize$y$center 
    beta <- beta/standardize$y$scale  
    gamma <- start$coefpath[,seq.int(length.out = q) + k]
    gamma <- t(apply(gamma, 1, standardize.coefficients, center = standardize$z$center, scale = standardize$z$scale))
    gamma[, grep("(Intercept)", colnames(gamma))] <- 
      gamma[, grep("(Intercept)", colnames(gamma))] - start$link$scale$linkfun(standardize$y$scale)  
    start$coefpath <- cbind(beta, gamma)
  } 
  

  if(is.character(dist)){
    ## distribution functions
    if(truncated) {
      ddist2 <- switch(dist, 
        "student"  = dtt, "gaussian" = dtnorm, "logistic" = dtlogis)
      sdist2 <- switch(dist, 
        "student"  = stt, "gaussian" = stnorm, "logistic" = stlogis)
    } else {
      ddist2 <- switch(dist, 
        "student"  = dct, "gaussian" = dcnorm, "logistic" = dclogis)
      sdist2 <- switch(dist, 
        "student"  = sct, "gaussian" = scnorm, "logistic" = sclogis)
    }
    ddist <- if(dist == "student") ddist2 else function(..., df) ddist2(...)
    sdist <- if(dist == "student") sdist2 else function(..., df) sdist2(...)


  } else { 
    ## for user defined distribution (requires list with ddist, sdist (optional)
    ## and hdist (optional), ddist, sdist, and hdist must be functions with
    ## arguments x, mean, sd, df, left, right, and log)
    ddist <- dist$ddist
    sdist <- if(is.null(dist$sdist)) NULL else  dist$sdist
    dist <- "user defined"
  }




 
  ## various fitted quantities (parameters, linear predictors, etc.)
  fitfun <- function(par, subset = 1:n) {
    beta <- par[seq.int(length.out = k)]
    gamma <- par[seq.int(length.out = q) + k]
    mu <- (drop(x %*% beta) + offset[[1L]])[subset]
    zgamma <- (drop(z %*% gamma) + offset[[2L]])[subset]
    sigma <- linkinv(zgamma)
    list(
      beta = beta,
      gamma = gamma,
      mu = mu,
      zgamma = zgamma,
      sigma = sigma
    )
  }

  ## objective function
  loglikfun <- function(par, subset = 1:n, sum = TRUE) {
    fit <- fitfun(par, subset)
    ll <- with(fit,  
        ddist(y[subset], mu, sigma, df = df, left = left[subset], right = right[subset], log = TRUE))
    if(sum) if(any(!is.finite(ll))) NaN else -sum(weights[subset] * ll)  
    else weights[subset] * ll
  }

  ## functions to evaluate gradient
  if(dfest | is.null(sdist)) {
    stop("score function required for boosting not available")
  } else { 
    gradfun <- function(par) {
      fit <- fitfun(par)
      grad <- with(fit, 
        sdist(y, mu, sigma, df = df, left = left, right = right))
      grad <- cbind(grad[,1], grad[,2] * mu.eta(fit$zgamma))
      return(-weights * grad)
    }
  }

  
  ## initialize intercepts
  if(is.null(start)) {
    par <- rep(0, k+q)
    coefpath <- par
    loglikpath <- -loglikfun(par)
    startit <- 1
  } else {
    par <- tail(start$coefpath, 1)
    coefpath <- start$coefpath
    loglikpath <- start$loglikpath
    startit <- start$iterations + 1 
  }


  ## actual boosting
  for(i in startit:maxit) {
    ## gradient of negative likelihood
    grad <- gradfun(par)
    ## location
    basefits <- .Call("mycov", -grad[,1],x, PACKAGE="crch")*n/sum(weights)
    par2 <- par
    minind <- which.max(abs(basefits))
    par2[minind] <- par2[minind] + nu*basefits[minind]
    ## scale
    basefits <- .Call("mycov", -grad[,2],z, PACKAGE="crch")*n/sum(weights)
    par3 <- par
    minind <- which.max(abs(basefits))
    par3[k + minind] <- par3[k + minind] + nu*basefits[minind]


    ## to compare mu and sigma improvements loglik must be used instead of rss
    ll3 <- -loglikfun(par3)
    ll2 <- -loglikfun(par2)
    par <- if(ll3 > ll2) par3 else par2
    coefpath <- rbind(coefpath, par)
    llnew <- max(ll3, ll2)
    loglikpath <- rbind(loglikpath, llnew)
  }

  ## cross validation to find optimum mstop
  if(mstop == "cv") {
    foldid <- control$foldid
    nfolds <- control$nfolds
    control$standardize <- FALSE
    control$mstop <- "max"
    control$start <- NULL
    if(is.null(foldid)) {
      foldid <- (sample(1:nfolds, size = length(y), replace = TRUE))
    } else {
      nfolds <- length(unique(foldid))
    }
  
    lltestall <- matrix(0, length(y), maxit + 1L)
    for(i in 1:nfolds) {
      train <- foldid != unique(foldid)[i]
      test <- foldid == unique(foldid)[i]
      boost <- crch.boost.fit(x = x[train, , drop = FALSE], y = y[train], 
        z = z[train, , drop = FALSE], left = left[train], right = right[train], 
        link.scale = link.scale, dist = dist, df = df, weights = weights[train], 
        offset = list(offset[[1L]][train], offset[[2L]][train]), 
        control = control, truncated = truncated)
      lltest <- apply(boost$coefpath, 1L, loglikfun, subset = test, sum = FALSE)
      lltest <- cbind(lltest, matrix(rep(lltest[,NCOL(lltest)], maxit + 1L - NCOL(lltest)), nrow = NROW(lltest)))
      lltestall[test,] <- lltest
    }

    mstopopt<- which.max(colMeans(lltestall))
  } else {
    mstopopt <- NULL
  }

      
  colnames(coefpath) <- c(colnames(x), paste("scale_", colnames(z), sep = ""))
  ## restandardize
  left  <-  left*standardize$y$scale + standardize$y$center
  right <- right*standardize$y$scale + standardize$y$center
  offset[[1L]] <- offset[[1L]]*standardize$y$scale
  beta <- coefpath[,seq.int(length.out = k), drop = FALSE]
  gamma <- coefpath[,seq.int(length.out = q) + k, drop = FALSE]
  beta[,] <- t(apply(beta, 1L, standardize.coefficients, 
    center = standardize$x$center, scale = standardize$x$scale, restandardize = TRUE))*standardize$y$scale
  beta[, grep("(Intercept)", colnames(beta))] <- 
    beta[, grep("(Intercept)", colnames(beta)), drop = FALSE] + standardize$y$center
  gamma[,] <- t(apply(gamma, 1, standardize.coefficients, 
    center = standardize$z$center, scale = standardize$z$scale, restandardize = TRUE))
  gamma[, grep("(Intercept)", colnames(gamma))] <- 
    gamma[, grep("(Intercept)", colnames(gamma)), drop = FALSE] + linkfun(standardize$y$scale)
  coefpath <- cbind(beta, gamma)
  y <- standardize.matrix(y, restandardize = TRUE)
  x <- standardize.matrix(x, restandardize = TRUE)
  z <- standardize.matrix(z, restandardize = TRUE)

  ## optimum stopping iterations
  mstopopt <- c("cv" = mstopopt, "max" = NROW(coefpath), 
      "aic" = which.min(2*rowSums(coefpath!=0) - 2* loglikpath),
      "bic" = which.min(log(nobs)*rowSums(coefpath!=0) - 2* loglikpath)) - 1

  ## coefficients, fitted values and likelihood for optimum stopping iteration
  par <- coefpath[mstopopt[[mstop]] + 1L,]
  fit <- fitfun(par)
  beta <- fit$beta
  gamma <- fit$gamma
  mu <- fit$mu
  sigma <- fit$sigma
  ll <- - loglikfun(par)

  names(beta) <- colnames(x)
  names(gamma) <- colnames(z)
  rownames(coefpath) <- c(0:(NROW(coefpath) - 1L))

  rval <- list(
    coefficients = list(location = beta, scale = gamma),
    df = df,
    residuals = y - mu,
    fitted.values = list(location = mu, scale = sigma),
    dist = dist,
    cens = list(left = left, right = right), 
    control = control,
    weights = if(identical(as.vector(weights), rep.int(1, n))) NULL else weights,
    offset = list(location = if(identical(offset[[1L]], rep.int(0, n))) NULL else 
      offset[[1L]],
      scale = if(identical(offset[[2L]], rep.int(0, n))) NULL else offset[[2L]]),
    n = n,
    nobs = nobs,
    loglik = ll,
    link = list(scale = linkobj),
    truncated = truncated,
#    converged = converged,
    coefpath = coefpath,
    loglikpath = loglikpath,
    iterations = NROW(coefpath) - 1L,
    stepsize = nu,
    mstop = mstop,
    mstopopt = mstopopt,
    standardize = standardize
  )
  
  class(rval) <- c("crch.boost", "crch")
  return(rval)
}


## Function to adjust stopping iteration of crch.boost objects
mstop.crch.boost <- function(object, mstop = NULL) {
  if(!is.null(mstop)) {
    ## adjust coefficients
    object$mstop <- mstop
    if(is.character(mstop)) {
      mstop <- match.arg(mstop, c("max", "aic", "bic", "cv"))
      mstop <- if(mstop == "max") object$iterations else object$mstopopt[mstop] 
    }
    k <- length(object$coefficient$location)
    q <- length(object$coefficient$scale)
    beta <- object$coefpath[mstop + 1,seq.int(length.out = k)]
    gamma <- object$coefpath[mstop + 1,seq.int(length.out = q) + k]
    names(gamma) <- names(object$coefficients$scale)
    object$coefficients <- list(location = beta, scale = gamma)
    ## remove or adjust some values
    object$residuals <- object$fitted.values$location <- object$fitted.values$scale <- NA
    object$loglik <- object$loglikpath[mstop + 1]
  }
  return(object)
}

coef.crch.boost <- function(object, model = c("full", "location", "scale", "df"), 
  mstop = NULL, zero.coefficients = FALSE, standardize = FALSE, ...) 
{
  model <- match.arg(model)
  object <- mstop.crch.boost(object, mstop = mstop)
  cf <- object$coefficients

  if(standardize){
    beta <- cf$location
    gamma <- cf$scale
    standardize <- object$standardize
    beta <- standardize.coefficients(beta, center = standardize$x$center, scale = standardize$x$scale)
    beta[grep("(Intercept)", colnames(beta))] <- beta[grep("(Intercept)", colnames(beta))] - standardize$y$center 
    beta <- beta/standardize$y$scale  
    gamma <- standardize.coefficients(gamma, center = standardize$z$center, scale = standardize$z$scale)
    gamma[grep("(Intercept)", colnames(gamma))] <- 
      gamma[grep("(Intercept)", colnames(gamma))] - object$link$scale$linkfun(standardize$y$scale)  
    cf <- list(location = beta, scale = gamma)
  }
  if(!zero.coefficients) {
    cf$location <- cf$location[cf$location != 0]
    cf$scale <- cf$scale[cf$scale != 0]
  }
  switch(model,
    "location" = {
      cf$location
    },
    "scale" = {
      cf$scale
    },
    "full" = {
      nam <- c(names(cf$location), paste("(scale)", names(cf$scale), sep = "_"))
      cf <- c(cf$location, cf$scale)
      names(cf) <- nam
      cf
    }
  )
}

print.crch.boost <- function(x, digits = max(3, getOption("digits") - 3), 
  mstop = NULL, zero.coefficients = FALSE, ...)
{
  x <- mstop.crch.boost(x, mstop = mstop)
  beta <- x$coefficients$location
  gamma <- x$coefficients$scale
  if(!zero.coefficients) {
    beta <- beta[beta != 0]
    gamma <- gamma[gamma != 0]
    prefix <- "Non-zero c"
  } else {
    prefix <- "C"
  }

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")
  cat(paste(prefix, "oefficients (location model):\n", sep = ""))
    print.default(format(beta, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
  cat(paste(prefix, "oefficients (scale model with ", x$link$scale$name, " link):\n", sep = ""))
    print.default(format(gamma, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")

  cat(paste("Distribution: ", x$dist, "\n", sep = ""))
  if(length(x$df)) {
    cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
  }
  cat("\n")

  invisible(x)
}

logLik.crch.boost <- function(object, mstop = NULL, ...) {
  object <- mstop.crch.boost(object)
  structure(object$loglik, df = sum(sapply(object$coefficients, function(x) sum(x != 0))), class = "logLik")
}


summary.crch.boost <- function(object, mstop = NULL, zero.coefficients = FALSE, ...)
{
  # mstop
  object <- mstop.crch.boost(object, mstop = mstop)
  mstop <- object$mstop   
  object$mstop <- if(mstop == "max") object$iterations else object$mstopopt[mstop] 
  # coefficients
  beta <- object$coefficients$location
  gamma <- object$coefficients$scale
  ## residuals
  object$residuals <- object$residuals/object$fitted.values$scale


  if(!zero.coefficients) {
    beta <- beta[beta != 0]
    gamma <- gamma[gamma != 0]
  }
  object$zero.coefficients <- zero.coefficients
  object$coefficients <- list(location = beta, scale = gamma)

  ## delete some slots
  object$fitted.values <- object$terms <- object$levels <- object$contrasts <- NULL

  ## return
  class(object) <- "summary.crch.boost"
  object
}

print.summary.crch.boost <- function(x, digits = max(3, getOption("digits") - 3),...) 
{
  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), sep = "\n")

  if(!all(is.na(x$residuals))) {
    cat(paste("\nStandardized residuals:\n", sep = ""))
    print(structure(round(as.vector(quantile(x$residuals)), digits = digits),
        .Names = c("Min", "1Q", "Median", "3Q", "Max")))
  }
  cat(paste("\nmaximum stopping iteration:", x$iterations, "\n"))
  cat(paste("\noptimum stopping iterations:\n", sep = ""))
  print(x$mstopopt)  
  cat("\n")
  prefix <- if(x$zero.coefficients) "C" else "Non-zero c"
  cat(paste(prefix, "oefficients after ", x$mstop, " boosting iterations:\n", sep = ""))
  cat(paste("Location model:\n", sep = ""))
  print.default(format(x$coefficients$location, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
   
  cat(paste("Scale model with", x$link$scale$name, "link:\n"))
  print.default(format(x$coefficients$scale, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
    
  cat(paste("Distribution: ", x$dist, "\n", sep = ""))
  if(length(x$df)) {
    cat(paste("Df: ", format(x$df, digits = digits), "\n", sep = ""))
  }
  cat("Log-likelihood:", formatC(x$loglik, digits = digits),
      "on", sum(sapply(x$coefficients, function(x) sum(x != 0))), "Df\n")
  cat("\n")
}

plot.crch.boost <- function(x, loglik = FALSE, 
  standardize = TRUE, which = c("both", "location", "scale"), mstop = NULL,
  coef.label = TRUE, col = NULL, ...) 
{
  which <- match.arg(which)
  k <- length(x$coefficient$location)
  q <- length(x$coefficient$scale)
  if(is.null(col)) col <- if(which == "both") c(1,2) else c(1,1)

  if(is.null(mstop)) {
    mstop <- x$mstop
  } else {
    mstop <- match.arg(mstop, c("max", "aic", "bic", "cv", "all", "no"))
  }
  if(mstop %in% c("no", "max")) {
    mstop <- NULL 
  } else {
    if(mstop == "all") mstop <- c("aic", "bic", "cv")
    mstop <- x$mstopopt[mstop] 
  }

  if(loglik) {
    coefupdate <- as.data.frame(apply(x$coefpath, 2, diff)!=0)
    coefupdate[["(Intercept)"]][rowSums(coefupdate) == 2] <- FALSE
    coefupdate[["scale_(Intercept)"]][rowSums(coefupdate) == 2] <- FALSE
    loglikpath <- x$loglikpath
    path <- rbind(0, diff(loglikpath)*coefupdate)
    path <- apply(path, 2, cumsum)
    colnames(path) <- c(names(x$coefficients$location), names(x$coefficients$scale))
    ylab <- "log-likelihood contribution"
  } else {
    beta <- x$coefpath[,seq.int(length.out = k), drop = FALSE]
    gamma <- x$coefpath[,seq.int(length.out = q) + k, drop = FALSE]
    if(standardize) {
      beta[,] <- t(apply(beta, 1, standardize.coefficients, 
        center = x$standardize$x$center, scale = x$standardize$x$scale))
      beta[, grep("(Intercept)", colnames(beta))] <- 
        beta[, grep("(Intercept)", colnames(beta))] - 
        x$standardize$y$center 
      beta <- beta/x$standardize$y$scale  
      gamma[,] <- t(apply(gamma, 1, standardize.coefficients, 
        center = x$standardize$z$center, scale = x$standardize$z$scale))
      gamma[, grep("(Intercept)", colnames(gamma))] <- 
        gamma[, grep("(Intercept)", colnames(gamma))] - 
        x$link$scale$linkfun(x$standardize$y$scale)  
      colnames(gamma) <- names(x$coefficients$scale)
    } 
    path <- cbind(beta, gamma)
    ylab <- if(standardize) "standardized coefficients" else "coefficients"
  }
  
  if(which == "location") {path <- path[,1:k]; q <- 0}
  if(which == "scale") {path <- path[,(k+1):(k+q)]; k <- 0}
  
  if(is.logical(coef.label)) {
    coef.label <- if(coef.label) {
        list(location = colnames(path)[1:k], scale = colnames(path)[(k+1):(q+k)])
    } else NULL
  }
  ## adjust figure margins
  rmar <- if(length(coef.label)) max(strwidth(colnames(x$coefpath), units = "in"))+0.5 else 0.42
  umar <- if(is.null(mstop)) 0.82 else 1.02
  par(mai = c(1.02, 0.82, umar, rmar))

  plot(ts(path, start = 0), plot.type = "single", ylab = ylab, xlab = "boosting iteration", 
    type = "s", col = rep(col, c(k,q)), ...)

  ## label paths
  if(length(coef.label)) {
    if(k!=0) axis(4, at = as.data.frame(tail(path, 1))[1:k][coef.label$location], 
      labels = coef.label$location, las = 1, col.axis = col[1], col.ticks = col[1])
    if(q!=0) axis(4, at = as.data.frame(tail(path, 1))[(k+1):(q+k)][coef.label$scale], 
      labels = coef.label$scale, las = 1, col.axis = tail(col, 1), col.ticks = tail(col,1))
  }

  ## vertical line at mstop
  if(!is.null(mstop)) {
    abline(v = mstop, lty = 2)
    axis(3, at = mstop, labels = names(mstop))
  }
}

predict.crch.boost <- function(object, newdata = NULL, mstop = NULL, 
  type = c("response", "location", "scale", "quantile"), na.action = na.pass, at = 0.5, 
  left = NULL, right = NULL, ...)
{
  object <- mstop.crch.boost(object, mstop = mstop)
  if(missing(newdata)) {
    if(!is.null(mstop)) stop("newdata has to be supplied for user defined mstop")
    predict.crch(object, type = type, na.action = na.action, at = at, left = left, right = right, ...)
  } else {
    predict.crch(object, newdata = newdata, type = type, na.action = na.action, 
      at = at, left = left, right = right, ...)
  }
}




  







