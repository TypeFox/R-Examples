################################################
## FRACTAL stationarity data constructor
## functions and corresponding methods
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class stationarity
## Constructor function: stationarity
## Methods:
##
##   as.list.stationarity
##   print.stationarity
##   print.summary.stationarity
##   summary.stationarity
##
################################################

###
# stationarity
###

"stationarity" <- function(x, n.taper=5,
  n.block=max(c(2, floor(logb(length(x), base=2)))), significance=0.05,
  center=TRUE, recenter=FALSE)
{
  # check input argument types and lengths
  checkScalarType(n.taper,"integer")
  checkScalarType(n.block,"integer")
  checkScalarType(significance,"numeric")
  checkScalarType(center,"logical")
  checkScalarType(recenter,"logical")

  series <- asVector(x)
  storage.mode(series) <- "double"

  if (is(x,"signalSeries"))
    series.name <- x@title
  else
    series.name <- deparse(substitute(x))

  if (!is.numeric(x) && !is(x, "signalSeries"))
    stop ("Input data must be a time series")

  # obtain sampling intercal of time series
  sampling.interval <- deltat(x)

  # perform argument checks
  n.sample <- length(x)

  if (n.taper < 5)
    stop("Must use at least 5 tapers in PSR stationarity tests")
  if (n.taper > n.sample)
    stop("Number of tapers is limited to the number of samples in the time series")
  if (n.sample <  2 * (12 * sampling.interval - 1))
    stop("The PSR stationarity test requires at least ",
      "2 * (12 * deltat(x) - 1) points in the time series")

  n.block.max <- ceiling(n.sample / (12.0 * sampling.interval - 1.0))

  if (n.block < 2 || n.block > n.block)
    stop(paste("n.block is out of range: 2 <= n.block <=", n.block.max))
  if (significance <= 0.0 || significance >= 1.0)
    stop("Significance out of range: 0 < significance < 1")

  z <- itCall( "RS_fractal_stationarity_priestley_subba_rao",
    series, as.double(sampling.interval),
    as.integer(n.taper), as.integer(n.block), as.double(significance),
    center, recenter)
    #COPY=rep(FALSE,7),
    #CLASSES=c("matrix","numeric","integer","integer","numeric","logical","logical"),
    #PACKAGE="ifultools")

  z[[3]] <- as.vector(z[[3]])

  names(z) <- c("test","anova","freq")

  dimnames(z[[1]]) <- list(c("Priestley-Subba Rao","Chi-square quantile"),
    c("T","I+R","T+I+R"))

  n.sample    <- length(series)
  n.frequency <- ncol(z[[2]]) - 1
  block.size  <- as.integer(floor(length(series) / n.block))
  deltaf      <- 1.0 / (block.size * sampling.interval)

  fourier.indices <- paste("f(", as.integer(z[[3]] / deltaf), ")", sep="")

  dimnames(z[[2]]) <- list(
    c(paste("block", seq(0, n.block - 1)), "col means"),
    c(fourier.indices, "row means"))

  names(z[[3]]) <- fourier.indices

  # form p-values

  dof <- c(n.block - 1, (n.block - 1) * (n.frequency - 1), (n.block - 1) * n.frequency)
  pvals <- 1 - pchisq(z[[1]][1,], dof)
  names(pvals) <- c("T","I+R","T+I+R")

  oldClass(z) <- "stationarity"

  attr(z, "method")            <- "Priestley-Subba Rao"
  attr(z, "pvals")             <- pvals
  attr(z, "n.sample")          <- n.sample
  attr(z, "sampling.interval") <- sampling.interval
  attr(z, "n.taper")           <- n.taper
  attr(z, "n.frequency")       <- n.frequency
  attr(z, "n.block")           <- n.block
  attr(z, "block.size")        <- block.size
  attr(z, "series.name")       <- series.name
  attr(z, "series")            <- series[seq(min(c(2048,length(series))))]
  attr(z, "significance")      <- significance
  attr(z, "deltaf")            <- deltaf
  attr(z, "fourier.indices")   <- fourier.indices
  attr(z, "stationarity")      <- z$test[1,] < z$test[2,]
  attr(z, "center")            <- center
  attr(z, "recenter")          <- recenter

  return(z)
}

###
# as.list.stationarity
###

"as.list.stationarity" <- function(x, ...)
  x[seq(along=x)]

###
# print.stationarity
###

"print.stationarity" <- function(x, justify="left", sep=":", digits=.Options$digits, ...)
{
  # define local functions
  "repchar" <- function(x, n) paste(rep(x, n), collapse="")

  xatt <- attributes(x)
  main <- paste(xatt$method, "stationarity Test for", xatt$series.name)

  z <- list(
    "Samples used"=xatt$n.sample,
    "Samples available"=xatt$block.size * xatt$n.block,
    "Sampling interval"=round(xatt$sampling.interval, digits),
    "SDF estimator"="Multitaper",
    "  Number of (sine) tapers"=xatt$n.taper,
    "  Centered"=xatt$center,
    "  Recentered"=xatt$recenter,
    "Number of blocks"=xatt$n.block,
    "Block size"=xatt$block.size,
    "Number of blocks"=xatt$n.block,
    "p-value for T"=xatt$pvals[1],
    "p-value for I+R"=xatt$pvals[2],
    "p-value for T+I+R"=xatt$pvals[3])

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  invisible(x)
}

###
# print.summary.stationarity
###

"print.summary.stationarity" <- function(x, ...)
{
  NextMethod("print")
  invisible(x)
}

###
# summary.stationarity
###

"summary.stationarity" <- function(object, ...)
{
  z <- as.list(object)
  oldClass(z) <- c("summary.stationarity","list")
  z
}
