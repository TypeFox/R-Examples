###################################################
## SAPA spectral density function constructor
## functions and corresponding methods
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class: SDF
## Constructor function: SDF
## Methods:
##
##   as.matrix.SDF
##   plot.SDF
##   print.SDF
##
## Other:
##
##   ACVS
##
###################################################

###
# SDF
###

"SDF" <- function(x, method="direct", taper.=NULL, window=NULL, n.taper=5, overlap=0.5,
  blocksize=NULL, single.sided=TRUE, sampling.interval=NULL, center=TRUE, recenter=FALSE,
  npad=2 * numRows(x))
{
	# check input arguments
  checkScalarType(method,"character")
  checkScalarType(n.taper,"integer")
  checkScalarType(overlap,"numeric")
  checkScalarType(single.sided,"logical")
  checkScalarType(center,"logical")
  checkScalarType(recenter,"logical")
  method <- match.arg(lowerCase(method),
    c("direct","lag window","wosa","multitaper"))

  if (is.list(x)){
    if (length(unique(unlist(lapply(x,numRows)))) > 1)
      stop("Inconsistent number of rows in input list objects")

    if (is.null(sampling.interval))
      sampling.interval <- deltat(x[[1]])

    series.name <- paste(names(x), collapse="|")
    if (nchar(series.name) == 0)
      series.name <- deparseText(substitute(x))

    x <- as.data.frame(x)
  }
  else if (isVectorAtomic(x)){

    if (is.null(sampling.interval))
      sampling.interval <- deltat(x)

    series.name <- ifelse1(is(x,"signalSeries"),
      x@title, deparseText(substitute(x)))

    x <- matrix(x, ncol=1)
  }
  else{
    series.name <- deparseText(substitute(x))
  }

  if (is.null(sampling.interval))
    sampling.interval <- deltat(x[,1])

  checkScalarType(npad,"integer")
  checkScalarType(sampling.interval,"numeric")

  if (sampling.interval <= 0)
    stop("Sampling interval must be positive")
  if (npad < 2)
    stop("npad must be at least 2")

  # define local functions
  "createTaper" <- function(taper., default.taper, n.sample, norm=TRUE, ...)
  {
    checkScalarType(default.taper,"character")
  	if (!is(taper.,"taper"))
      taper. <- taper(type=ifelse1(is.null(taper.), default.taper, taper.),
        n.sample=n.sample, normalize=norm, ...)

    # check taper attributes
    if (attr(taper., "n.sample") != n.sample)
      stop("Taper vector must contain ", n.sample, " elements")

    taper.
  }

  # convert to a double matrix
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n.sample <- numRows(x)

  if ((method != "lag window" && npad < n.sample) ||
      (method == "lag window" && npad < 2 * n.sample)){

    M <- floor(npad/2) + 1
    decimation.factor <- ceiling((n.sample + 1)/M)
    npad <- 2 * (decimation.factor * M) - 2
  }
  else
    decimation.factor <- 1

  if (method == "direct"){

    # issue periodogram by default
    taper. <- createTaper(taper., "rectangle", n.sample)

    S <- itCall("RS_fractal_spectral_density_function_direct",
      x, as.vector(taper.),
      as.logical(center), as.logical(recenter), TRUE, as.integer(npad))
        #,
      #COPY=rep(FALSE,6),
      #CLASSES = c(rep("matrix", 2), rep("logical",3), "integer"),
      #PACKAGE="ifultools")
  }
  else if (method == "lag window"){

    # taper. is the lag window and window is the taper used
    # in the interim direct SDF estimation
    taper. <- createTaper(taper., "parzen", n.sample, cutoff = n.sample %/% 2)
    window <- createTaper(window, "hanning", n.sample)

    S <- itCall("RS_fractal_spectral_density_function_lag_window",
      x, as.vector(taper.), as.vector(window),
      as.logical(center), as.logical(recenter), TRUE, as.integer(npad))
      #,
      #COPY=rep(FALSE,7),
      #CLASSES = c(rep("matrix", 3), rep("logical",3), "integer"),
      #PACKAGE="ifultools")
  }
  else if (method == "wosa"){

    if (is.null(blocksize))
      blocksize <- as.integer(n.sample %/% 4)
    checkScalarType(blocksize,"integer")

    taper. <- createTaper(taper., "hanning", blocksize)

    if (overlap < 0.0 || overlap >= 1)
      stop("Overlap fraction must be on 0 <= overlap < 1")

    S <- itCall("RS_fractal_spectral_density_function_wosa",
      x, as.vector(taper.), as.numeric(overlap),
      as.logical(center), as.logical(recenter), TRUE, as.integer(npad))
        #,
      #COPY=rep(FALSE,7),
      #CLASSES = c(rep("matrix", 2), "numeric", rep("logical",3), "integer"),
      #PACKAGE="ifultools")
  }
  else if (method == "multitaper"){

    taper. <- createTaper(taper., "sine", n.sample, n.taper=min(c(n.sample, n.taper)))

    S <- itCall("RS_fractal_spectral_density_function_multitaper",
      x, as.matrix(taper.),
      as.logical(center), as.logical(recenter), TRUE, as.integer(npad))
  #,
      #COPY=rep(FALSE,6),
      #CLASSES = c(rep("matrix", 2), rep("logical", 3), "integer"),
      #PACKAGE="ifultools")
  }
  else
    stop("SDF estimator is unsupported")

  # convert to a vector and adjust for sampling interval
  S <- S * sampling.interval

  # define frequency resolution
  M  <- npad %/% 2 + 1
  df <- 1 / npad / sampling.interval
  f  <- seq(0, M-1) * df

  # decimate results if necessary
  if (decimation.factor > 1){
	ix <- seq(1, M, by=decimation.factor)
    S  <- S[ix,]
    f  <- f[ix]
  }

  # adapt for double-sided if requested
  if (!single.sided){
	ix <- seq(along=f)
	ix <- c(rev(ix[-1]),ix)
    S  <- S[ix,]
    f  <- c(-rev(f[-1]),f)
  }

  # convert to all real values if number of series is unity
  if (numCols(x) == 1)
    S <- as.vector(Re(S))

  # create cross-SDF labels
  n.series <- numCols(x)
  S.labels <- NULL
   for (i in seq(n.series))
     for (j in seq(i,n.series))
        S.labels <- c(S.labels, paste("S",i,j,sep=""))

  # assign class and attributes
  oldClass(S) <- "SDF"

  attr(S, "frequency")      <- f
  attr(S, "method")         <- method
  attr(S, "taper")          <- taper.
  attr(S, "lag.window")     <- window
  attr(S, "n.taper")        <- n.taper
  attr(S, "wosa.blocksize") <- blocksize
  attr(S, "wosa.overlap")   <- overlap
  attr(S, "wosa.nblock")    <- 1 + floor((n.sample - blocksize) / (blocksize * (1 - overlap)))
  attr(S, "deltat")         <- sampling.interval
  attr(S, "deltaf")         <- abs(diff(f[1:2]))
  attr(S, "single.sided")   <- single.sided
  attr(S, "series.name")    <- series.name
  attr(S, "n.sample")       <- n.sample
  attr(S, "center")         <- center
  attr(S, "recenter")       <- recenter
  attr(S, "n.series")       <- n.series
  attr(S, "labels")         <- S.labels
  attr(S, "n.sdf")          <- length(S.labels)
  S
}

###
# as.matrix.SDF
###

"as.matrix.SDF" <- function(x, ...)
   matrix(as.vector(x), ncol=attr(x,"n.sdf"), ...)

###
# plot.SDF
###

"plot.SDF" <- function(x, xscale="linear", yscale="dB", type="l",
  xlab="FREQUENCY (Hz)", ylab=NULL, plot.mean=!attr(x,"center"), n.plot=3, FUN=Mod, add=FALSE,
  col="black", ...)
{
  checkScalarType(type,"character")
  checkScalarType(xlab,"character")
  checkScalarType(plot.mean,"logical")
  checkScalarType(n.plot,"integer")
  checkScalarType(add,"logical")
  if (!is(FUN,"function"))
    stop("FUN must be a function")

  xatt <- attributes(x)
  S    <- as.matrix(x)
  f    <- xatt$frequency

  if (!plot.mean){
    nix <- ifelse1(xatt$single.sided, 1, ceiling(numRows(x)/2))
    f   <- f[-nix]
    S   <- S[-nix,,drop=FALSE]
  }

  S.old <- S

  # ensure a supported FUN if data is complex
  FUN.str <- match.arg(deparse(substitute(FUN)), c("Mod","Im","Re","Arg"))
  if (is.complex(S)){
    if (is.null(ylab))
      ylab <- paste(properCase(xatt$method))
  }
  else{
    if (is.element(FUN.str, c("Im","Arg")))
      warning("SDF is purely real: coercing FUN to Mod")
    FUN <- Mod
    FUN.str <- "Mod"
    if (is.null(ylab))
      ylab <- paste(properCase(xatt$method),"SDF")
  }

  # convert data
  S.fun <- FUN(S)

  # check FUN and scaling
  if (any(S.fun < 0) && yscale != "linear"){
   yscale <- "linear"
   warning("Obtained negative values after applying specified FUN function to SDF. Coercing yscale to linear")
  }

  S     <- scaleData(S.fun, yscale)
  f     <- as.vector(scaleData(f, xscale))
  xlab  <- paste(attr(f,"scalestr"), xlab)

  if (all(is.na(S))){
   FUN.str <- "Re"
   S <- Re(S.old)
   warning("Obtained all NA values after applying specified FUN function to SDF. Using FUN=Re instead")
  }

  # specify plot grid
  if (!add){
    old.mfrow <- par(mfrow=c(1,1))
    on.exit(par(old.mfrow))
  }

  if (xatt$n.series == 1){

    ylab <- paste(attr(S,"scalestr"), ylab)
    plot(x=f, y=S, type=type, xlab=xlab, ylab=ylab, col=col)
  }
  else{

   S.labels <- paste(attr(S,"scalestr"), FUN.str, xatt$labels)
   iz <- seq(n.plot)
   main <- list(text=paste(ylab, "Cross SDF Estimates"), adj=0, cex=0.7)

   for (i in seq(numCols(S) %/% n.plot)){

      Z <- as.data.frame(S[,iz])

      stackPlot(f, Z, xlab=xlab, ylab=S.labels[iz], add=add, main=main, col=col)
      iz <- iz + n.plot
   }

   ileft <- numCols(S) %% n.plot
   if (ileft){
     iz <- numCols(S) - rev(seq(ileft))
     Z  <- as.data.frame(S[,iz])

     stackPlot(f, Z, xlab=xlab, ylab=S.labels[iz], add=add, main=main, col=col)
   }
 }

 invisible(NULL)
}

###
# print.SDF
###

"print.SDF" <- function(x, justify="left", sep=":", ...)
{
  xatt   <- attributes(x)
  taper. <- xatt$taper
  sub    <- xatt$lag.window
  cross  <- xatt$n.series > 1
  wosa   <- xatt$method == "wosa"
  mult   <- xatt$method == "multitaper"
  main   <- paste(ifelse1(cross,"Cross-", ""),
    "Spectral Density Function estimation for ", xatt$series.name, sep="")

  z <- list(
    "Cross-SDF labels"=ifelse1(cross, xatt$labels, NULL),
    "Length of series"=xatt$n.sample,
    "Sampling interval"=xatt$deltat,
    "Frequency resolution (Hz)"=xatt$deltaf,
    "Centered"=xatt$center,
    "Recentered"=xatt$recenter,
    "Single-sided"=xatt$single.sided,
    "Method"=properCase(xatt$method),
    "  Block size"=ifelse1(wosa, xatt$wosa.blocksize, NULL),
    "  Overlap"=ifelse1(wosa, xatt$wosa.overlap, NULL),
    "  Number of blocks"=ifelse1(wosa, xatt$wosa.nblock, NULL),
    "Number of tapers"=ifelse1(mult, xatt$n.taper, NULL)
   )

  prettyPrintList(z, header=main, sep=sep, justify=justify, ...)

  if (is(taper.,"taper"))
    print(taper., pre="  ")
  if (is(sub,"taper"))
    print(sub, pre="  ")

  invisible(x)
}

###
# ACVS
###

"ACVS" <- function(x, biased=TRUE, center=TRUE)
{
  as.vector(itCall("RS_math_acvs",
    as.numeric(x), as.logical(biased), as.logical(center))
    #,
    #COPY = rep(FALSE,3),
    #CLASSES = c("numeric",rep("logical",2)),
    #PACKAGE="ifultools")
)
}
