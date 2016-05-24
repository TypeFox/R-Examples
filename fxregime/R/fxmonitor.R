fxmonitor <- function(formula, data, start, end = 3, alpha = 0.05, meat. = NULL)
{
  ## formula processing
  if(missing(formula)) formula <- colnames(data)[1]
  if(is.character(formula)) formula <- as.formula(paste(formula,
    "~", paste(colnames(data)[colnames(data) != formula], collapse = " + ")))

  ## historical data
  if(is.character(start)) start <- as.Date(start)
  monitor <- end(window(data, end = start-1))
  n <- NROW(window(data, end = monitor))
  hist.fm <- fxlm(formula, data = window(data, end = monitor))
  sigma2 <- coef(hist.fm)["(Variance)"]

  ## full data
  mf <- model.frame(formula, data = as.data.frame(data))
  y <- model.response(mf)
  X <- model.matrix(formula, data = as.data.frame(data))
  attr(X, "assign") <- NULL

  ## estimating functions
  y.pred <- predict(hist.fm, as.data.frame(X))
  mscore <- as.vector(y - y.pred) * X
  mscore <- cbind(mscore, ((y - y.pred)^2 - sigma2))

  ## cumulative estimating functions
  mscore <- mscore/sqrt(n)
  hist.score <- mscore[1:n,]
  J12 <- if(!is.null(meat.)) meat.(hist.fm) else crossprod(hist.score)
  J12 <- root.matrix(J12)

  mscore <- apply(as.matrix(mscore), 2, cumsum)
  mscore <- zoo(t(chol2inv(chol(J12)) %*% t(mscore)), order.by = index(data))
  colnames(mscore) <- names(coef(hist.fm))

  ## critical value table (Table III from Zeileis et al, JAE)
  cv <- matrix(c(
    1.159, 1.329, 1.430, 1.472, 1.502, 1.541, 1.567,
    1.253, 1.445, 1.544, 1.589, 1.619, 1.668, 1.688,
    1.383, 1.590, 1.695, 1.753, 1.789, 1.838, 1.860,
    1.467, 1.688, 1.793, 1.861, 1.899, 1.961, 1.964,
    1.568, 1.814, 1.939, 2.006, 2.046, 2.090, 2.128,
    1.616, 1.896, 2.022, 2.076, 2.131, 2.159, 2.219,
    1.680, 1.997, 2.103, 2.177, 2.226, 2.257, 2.311,
    1.801, 2.114, 2.217, 2.301, 2.397, 2.380, 2.454,
    1.976, 2.300, 2.423, 2.525, 2.573, 2.597, 2.650,
    2.118, 2.478, 2.599, 2.712, 2.812, 2.766, 2.888,
    2.435, 2.789, 2.973, 3.288, 3.226, 3.230, 3.401), nrow = 7)
  cv_row <- c(2:6, 8, 10)
  cv_col <- c(0.2, 0.15, 0.1, 0.075, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001)

  ## choose Bonferroni-corrected critical value
  alpha <- 1 - (1 - alpha)^(1/ncol(mscore))
  critval <- apply(cv, 2, function(x) approx(cv_row, x, end, rule = 2)$y)
  critval <- approx(cv_col, critval, alpha, rule = 2)$y
  
  rval <- list(process = mscore,
               n = n,
	       formula = formula,
	       data = data,
	       monitor = monitor,
	       critval = critval,
	       J12 = J12)
  class(rval) <- "fxmonitor"
  return(rval)
}

plot.fxmonitor <- function(x,
  which = NULL, aggregate = NULL,
  ylim = NULL, xlab = "Time", ylab = "Empirical fluctuation process",
  main = "Monitoring of FX model", ...)
{
  n <- x$n
  proc <- x$process
  if(is.null(which)) which <- 1:NCOL(proc)
  if(is.null(aggregate)) aggregate <- length(which) > 1
  
  proc <- proc[,which]
  dummy <- zoo((n+1):NROW(proc)/n, index(proc)[-(1:n)]) * x$critval

  if(aggregate) {
    maxabs <- zoo(apply(abs(coredata(proc)), 1, max), order.by = index(proc))
    if(is.null(ylim)) ylim <- c(0, max(c(max(maxabs), max(dummy))))
    plot(maxabs, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
    abline(h = 0)
    abline(v = x$monitor, lty = 2)
    lines(dummy, col = 2)
  } else {
    if(missing(ylim)) ylim <- c(min(min(proc), min(-dummy)), max(max(proc), max(dummy)))
    if(missing(ylab) && NCOL(proc) > 1) ylab <- colnames(proc)
    mypanel <- function(z, ...) {
      lines(z, ...)
      abline(h = 0)
      abline(v = x$monitor, lty = 2)
      lines(-dummy, col = 2)
      lines(dummy, col = 2)
    }
    plot(proc, panel = mypanel, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  }
  invisible(x)
}

print.fxmonitor <- function(x, ...)
{
  n <- x$n
  proc <- x$process
  maxabs <- apply(abs(coredata(proc)), 1, max)[-(1:n)]
  dummy <- ((n+1):NROW(proc)/n) * x$critval
  breakpoint <- if(any(maxabs > dummy)) min(which(maxabs > dummy)) + n else NA
  form <- as.character(x$formula)
  form <- paste(form[2], form[1], form[3])

  cat("Monitoring of FX model\n\n")
  cat(paste("Formula:", form, "\n"))
  cat(paste("History period:", start(proc), "to", x$monitor, "\n"))
  cat(paste("Break detected:", ifelse(is.na(breakpoint), "none", as.character(index(proc)[breakpoint])), "\n"))
  
  invisible(x)
}

breakpoints.fxmonitor <- function(obj, ...)
{
  n <- obj$n
  proc <- obj$process
  maxabs <- apply(abs(coredata(proc)), 1, max)[-(1:n)]
  dummy <- ((n+1):NROW(proc)/n) * obj$critval
  breakpoint <- if(any(maxabs > dummy)) min(which(maxabs > dummy)) + n else NA
  return(breakpoint)
}

breakdates.fxmonitor <- function(obj, ...)
{
  rval <- breakpoints(obj)
  if(!is.na(rval)) rval <- index(obj$process)[rval]
  return(rval)
}
