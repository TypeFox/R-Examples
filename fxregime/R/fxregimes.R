## partition FX models based on normal loglik
## (wrapper to the currently unexported gbreakpoints)
fxregimes <- function(formula, data, ..., hpc = c("none", "foreach")) {
  if(missing(formula)) formula <- colnames(data)[1]
  if(is.character(formula)) formula <- as.formula(paste(formula,
    "~", paste(colnames(data)[colnames(data) != formula], collapse = " + ")))
  rval <- gbreakpoints(formula, data = data, order.by = time(data), hpc = hpc, ...)
  rval$formula <- formula
  rval$data <- data
  class(rval) <- c("fxregimes", class(rval))
  return(rval)
}

print.fxregimes <- function(x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nFX model: ", deparse(x$formula), "\n\n",
        paste("Minimum", x$icname, "partition\n"),
	"Number of regimes: ", length(x$breakpoints[!is.na(x$breakpoints)]) + 1, "\n",
	"Breakdates: ", paste(breakdates(x), collapse = ", "), "\n\n",
	sep = "")
    invisible(x)
}

index.fxregimes <- time.fxregimes <- function(x, ...) index(x$data)

lines.fxregimes <- function(x, breaks = NULL, lty = 2, ...)
  abline(v = breakdates(x, breaks = breaks), lty = lty, ...)


## obtain FX models for each segment
## (plus standard extractor functions)
refit.fxregimes <- function(object, breaks = NULL, ...) {
  if(is.null(breaks)) breaks <- time(object$data)[object$breakpoints]
    else if(length(breaks) == 1 & is.numeric(breaks) & !inherits(breaks, "Date"))
      breaks <- breakdates(object, breaks = breaks)
  breaks <- breaks[!is.na(breaks)]
  
  if(length(breaks) < 1) {
    sbp <- start(object$data)
    ebp <- end(object$data)
  } else {  
    sbp <- c(start(object$data), sapply(breaks, function(z)
      index(object$data)[min(which(index(object$data) > z))]))
    ebp <- c(breaks, end(object$data))
  }
  
  rval <- lapply(1:length(sbp), function(i)
    fxlm(object$formula, data = window(object$data, start = sbp[i], end = ebp[i]), ...))
  names(rval) <- paste(as.character(sbp), as.character(ebp), sep = "--")
  return(rval)  
}

coef.fxregimes <- function(object, breaks = NULL, ...)
  t(sapply(refit(object, breaks = breaks, ...), coef))

fitted.fxregimes <- function(object, breaks = NULL, ...) {
  rval <- as.vector(unlist(lapply(refit(object, breaks = breaks, ...), fitted)))
  rval <- zoo(rval, index(object))
  return(rval)
}

residuals.fxregimes <- function(object, breaks = NULL, ...) {
  rval <- as.vector(unlist(lapply(refit(object, breaks = breaks, ...), residuals)))
  rval <- zoo(rval, index(object))
  return(rval)
}

confint.fxregimes <- function(object, parm = NULL, level = 0.95, breaks = NULL,  meat. = NULL, ...)
{
  ## parameters: level, breaks
  a2 <- (1 - level)/2
  if(!is.null(parm) & !is.null(breaks))
    warning("`parm' and `breaks' are both specified: `breaks' is used")
  else
    if(!is.null(parm)) breaks <- parm

  ## extract estimates
  bp <- breakpoints(object, breaks = breaks)$breakpoints
  if(any(is.na(bp))) stop("cannot compute confidence interval when `breaks = 0'")
  
  ## set up intervals
  nbp <- length(bp)
  upper <- rep(0, nbp)
  lower <- rep(0, nbp)
  bp <- c(0, bp, object$nobs)

  ## fitted models
  rf <- refit(object, breaks = breaks)
  cf <- lapply(rf, coef)
  Q <- lapply(rf, function(x) solve(bread(x)))
  Omega <- if(is.null(meat.)) Q else lapply(rf, meat.)

  ## utility functions
  myfun <- function(x, level = 0.975, xi = 1, phi1 = 1, phi2 = 1)
    (pargmaxV(x, xi = xi, phi1 = phi1, phi2 = phi2) - level)
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
 
    p0 <- pargmaxV(0, phi1 = phi1, phi2 = phi2, xi = xi)
    if(is.nan(p0) || p0 < a2 || p0 > (1-a2)) {
      warning(paste("Confidence interval", as.integer(i),
        "cannot be computed: P(argmax V <= 0) =", round(p0, digits = 4)))
      upper[i-1] <- NA
      lower[i-1] <- NA
    } else {
      ub <- lb <- 0
      while(pargmaxV(ub, phi1 = phi1, phi2 = phi2, xi = xi) < (1 - a2)) ub <- ub + 1000
      while(pargmaxV(lb, phi1 = phi1, phi2 = phi2, xi = xi) > a2) lb <- lb - 1000

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
  class(RVAL) <- "confint.fxregimes"
  return(RVAL)
}

breakdates.confint.fxregimes <- function(obj, format.times = FALSE, ...)
{
  data.frame(lapply(as.data.frame(obj$confint), function(x) obj$index[x]), check.names = FALSE)
}

print.confint.fxregimes <- function(x, ...)
{
  nbp <- nrow(x$confint)
  cat("\n\t Confidence intervals for breakpoints")
  cat(paste("\n\t of optimal ", nbp + 1, "-segment partition: \n\n", sep = ""))
  cat("Call:\n")
  print(x$call)
  cat("\nBreakpoints at observation number:\n")
  print(x$confint, quote = FALSE)
  cat("\nCorresponding to breakdates:\n")
  print(breakdates(x), quote = FALSE)
}

lines.confint.fxregimes <- function(x, col = 2, angle = 90, length = 0.05,
  code = 3, at = NULL, breakpoints = TRUE, ...)
{
  nbp <- nrow(x$confint)
  x <- breakdates(x)
  if(breakpoints) abline(v = x[,2], lty = 2)
  if(is.null(at)) {
    at <- par("usr")[3:4]
    at <- diff(at)/1.08 * 0.02 + at[1]
  }
  if(length(at) < nbp) at <- rep(at, length.out = nbp)
  arrows(x[,1], at, x[,3], at, col = col, angle = angle, length = length, code = code, ...)
}
