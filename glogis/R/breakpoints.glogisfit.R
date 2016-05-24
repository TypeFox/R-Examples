###########################
## Breakpoint estimation ##
###########################

## This is somewhat experimental...
## glue code for gbreakpoints() in fxregime
breakpoints.glogisfit <- function(obj, h = 0.15, breaks = NULL, ic = "LWZ", hpc = "none", ...)
{
  stopifnot(require("fxregime"))
  dat <- data.frame(x = as.vector(obj$x))
  glogisfit0 <- function(formula, data, ...) glogisfit.default(data$x, fixed = obj$fixed, hessian = FALSE, ...)
  rval <- fxregime:::gbreakpoints(x ~ 1, data = dat, order.by = time(obj$x), fit = glogisfit0,
    h = h, breaks = breaks, ic = ic, hpc = hpc)
  rval$null <- obj
  class(rval) <- c("breakpoints.glogisfit", class(rval))
  return(rval)
}

refit.breakpoints.glogisfit <- function(object, ...) {
  bf <- strucchange::breakfactor(object, ...)
  rval <- tapply(object$null$x, bf, glogisfit, fixed = object$null$fixed)
  names(rval) <- paste(tapply(format(object$index), bf, head, 1), "--",
    tapply(format(object$index), bf, tail, 1), sep = "")
  return(rval)
}

coef.breakpoints.glogisfit <- function(object, log = TRUE, ...) {
  rf <- fxregime::refit(object, ...)
  t(sapply(rf, coef, log = log))
}

fitted.breakpoints.glogisfit <- function(object,
  type = c("mean", "variance", "skewness"), ...)
{
  type <- as.vector(sapply(type, match.arg, choices = c("mean", "variance", "skewness")))
  rf <- fxregime::refit(object, ...)
  mom <- t(sapply(rf, "[[", "moments"))
  rval <- mom[strucchange::breakfactor(object, ...), type]
  if(inherits(object$null$x, "zoo")) rval <- zoo(rval, time(object$null$x))
  if(inherits(object$null$x, "ts")) rval <- ts(rval, start = start(object$null$x), frequency = frequency(object$null$x))
  return(rval)
}

index.breakpoints.glogisfit <- function(x, ...) x$index

confint.breakpoints.glogisfit <- function(object, parm = NULL, level = 0.95, breaks = NULL,  meat. = NULL, ...)
{
  ## parameters: level, breaks
  a2 <- (1 - level)/2
  if(!is.null(parm) & !is.null(breaks))
    warning("`parm' and `breaks' are both specified: `breaks' is used")
  else
    if(!is.null(parm)) breaks <- parm

  ## extract estimates
  bp <- strucchange::breakpoints(object, breaks = breaks)$breakpoints
  if(any(is.na(bp))) stop("cannot compute confidence interval when `breaks = 0'")
  
  ## set up intervals
  nbp <- length(bp)
  upper <- rep(0, nbp)
  lower <- rep(0, nbp)
  bp <- c(0, bp, object$nobs)

  ## fitted models
  rf <- fxregime::refit(object, breaks = breaks)
  cf <- lapply(rf, coef)
  Q <- lapply(rf, function(x) solve(bread(x)))
  Omega <- if(is.null(meat.)) Q else lapply(rf, meat.)

  ## utility functions
  myfun <- function(x, level = 0.975, xi = 1, phi1 = 1, phi2 = 1)
    (strucchange::pargmaxV(x, xi = xi, phi1 = phi1, phi2 = phi2) - level)
  myprod <- function(delta, mat) as.vector(crossprod(delta, mat) %*% delta)

  ## loop over breaks
  for(i in 1:nbp)
  {
    delta <- cf[[i+1]] - cf[[i]]            
    Oprod1 <- myprod(delta, Omega[[i]])
    Oprod2 <- myprod(delta, Omega[[i+1]])
    Qprod1 <- myprod(delta, Q[[i]])
    Qprod2 <- myprod(delta, Q[[i+1]])

    xi <- Qprod2/Qprod1
    phi1 <- sqrt(Oprod1/Qprod1)
    phi2 <- sqrt(Oprod2/Qprod2)
 
    p0 <- strucchange::pargmaxV(0, phi1 = phi1, phi2 = phi2, xi = xi)
    if(is.nan(p0) || p0 < a2 || p0 > (1-a2)) {
      warning(paste("Confidence interval", as.integer(i),
        "cannot be computed: P(argmax V <= 0) =", round(p0, digits = 4)))
      upper[i-1] <- NA
      lower[i-1] <- NA
    } else {
      ub <- lb <- 0
      while(strucchange::pargmaxV(ub, phi1 = phi1, phi2 = phi2, xi = xi) < (1 - a2)) ub <- ub + 1000
      while(strucchange::pargmaxV(lb, phi1 = phi1, phi2 = phi2, xi = xi) > a2) lb <- lb - 1000

      upper[i] <- uniroot(myfun, c(0, ub), level = (1-a2), xi = xi, phi1 = phi1, phi2 = phi2)$root
      lower[i] <- uniroot(myfun, c(lb, 0), level = a2, xi = xi, phi1 = phi1, phi2 = phi2)$root
    
      upper[i] <- upper[i] * phi1^2 / Qprod1
      lower[i] <- lower[i] * phi1^2 / Qprod1
    }
  }
  
  ## collect results
  bp <- bp[-c(1, nbp+2)]
  bp <- cbind(bp - ceiling(upper), bp, bp - floor(lower))
  a2 <- round(a2 * 100, digits = 1)
  colnames(bp) <- c(paste(a2, "%"), "breakpoints", paste(100 - a2, "%"))
  rownames(bp) <- 1:nbp
  RVAL <- list(confint = bp,
               nobs = object$nobs,
	       npar = object$npar,
	       call = match.call(),
               index = index(object))
  class(RVAL) <- "confint.breakpoints.glogisfit"
  return(RVAL)
}

breakdates.confint.breakpoints.glogisfit <- function(obj, format.times = FALSE, ...)
{
  data.frame(lapply(as.data.frame(obj$confint), function(x) obj$index[x]), check.names = FALSE)
}

print.confint.breakpoints.glogisfit <- function(x, ...)
{
  nbp <- nrow(x$confint)
  cat("\n\t Confidence intervals for breakpoints")
  cat(paste("\n\t of optimal ", nbp + 1, "-segment partition: \n\n", sep = ""))
  cat("Call:\n")
  print(x$call)
  cat("\nBreakpoints at observation number:\n")
  print(x$confint, quote = FALSE)
  cat("\nCorresponding to breakdates:\n")
  print(strucchange::breakdates(x), quote = FALSE)
}

lines.confint.breakpoints.glogisfit <- function(x, col = 2, angle = 90, length = 0.05,
  code = 3, at = NULL, breakpoints = TRUE, ...)
{
  nbp <- nrow(x$confint)
  x <- strucchange::breakdates(x)
  if(breakpoints) abline(v = x[,2], lty = 2)
  if(is.null(at)) {
    at <- par("usr")[3:4]
    at <- diff(at)/1.08 * 0.02 + at[1]
  }
  if(length(at) < nbp) at <- rep(at, length.out = nbp)
  arrows(x[,1], at, x[,3], at, col = col, angle = angle, length = length, code = code, ...)
}
