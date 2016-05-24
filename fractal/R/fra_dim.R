################################################
## FRACTAL nonlinear dynamics functionality
##
##    corrDim (d2)
##    infoDim (d1)
##    embedSeries
##    FNN
##    FNS
##    lyapunov
##    poincareMap
##    spaceTime
##    timeLag
##
################################################

###
# corrDim
###

"corrDim" <- function(x, dimension=5,
  tlag=timeLag(x, method="acfdecor"), olag=0, resolution=2)
{
  # check input argument types and lengths
  checkScalarType(dimension,"integer")
  checkScalarType(olag,"integer")
  checkScalarType(resolution,"numeric")

  # obtain series name
  if (is(x,"embedSeries"))
    series.name <- attr(x,"series.name")
  else
    series.name <- deparseText(substitute(x))

  # check to see if input is an embedding or a time series
  if (!is.numeric(x) && !is(x,"signalSeries"))
    stop ("Input data must be a numeric vector or an object of class \"signalSeries\"")

  # embed the data if not already embedded
  if (!is(x,"embedSeries"))
    x <- embedSeries(x, dimension=dimension, tlag=tlag, series.name=series.name)

  # obtain attributes
  xatt      <- attributes(x)
  dimension <- max(xatt$emb.dim)
  n.embed   <- xatt$n.embed
  tlag      <- xatt$tlag

  # calculate the correlation summation
  z <- itCall( "RS_fractal_dimension_correlation_summation",
    as.matrix(x), dimension, tlag, olag, resolution)
    #COPY=rep(FALSE,5), #CLASSES=c("matrix",rep("integer",4)),
    #PACKAGE="ifultools")

  # correlation summations greater than unity
  # are a flag for scales that contain no data.
  # replace these with NAs
  z[[1]][z[[1]] > 1] <- NA

  if (length(z[[2]]) < 5)
    stop("Insufficient number of scales calculated to form C2 statistics")

  chaoticInvariant(z,
    dimension   = seq(dimension),
    n.embed     = n.embed,
    n.reference = n.embed,
    n.neighbor  = NA,
    tlag        = tlag,
    olag        = olag,
    resolution  = resolution,
    series.name = series.name,
    series      = as.vector(x[,1]),
    ylab        = "log2(C2)",
    xlab        = "log2(scale)",
    metric      = Inf,
    invariant   = "correlation dimension")
}

###
# infoDim
###

"infoDim" <- function(x, dimension=5, tlag=NULL,
  olag=0, n.density=100, metric=Inf,
  max.neighbors=as.integer(min(c(round(length(x) / 3), 100))),
  n.reference=as.integer(round(length(x) / 20)))
{
  if (!is.numeric(x) && !is(x,"signalSeries"))
    stop ("Input data must be numeric or an object of class")
  if (is.null(tlag))
    tlag <- timeLag(x, method="acfdecor")

  # check input argument types and lengths
  checkScalarType(dimension,"integer")
  checkScalarType(metric,"numeric")
  checkScalarType(olag,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(n.reference,"integer")
  checkScalarType(n.density,"integer")
  checkScalarType(max.neighbors,"integer")

  # obtain series name
  series.name <- deparseText(substitute(x))

  if (dimension < 2)
    stop("Maximal embedding dimensio must be greater than two")


  z <- itCall( "RS_fractal_dimension_information",
    as.numeric(x), as.integer(dimension), as.integer(tlag), as.integer(olag), as.integer(n.density),
    mutilsDistanceMetric(metric), as.integer(max.neighbors), as.integer(n.reference))
    #COPY=rep(FALSE,8), #CLASSES=c("matrix", rep("integer", 7)),
    #PACKAGE="ifultools")

  chaoticInvariant(z,
    dimension   = seq(2,dimension),
    n.embed     = as.integer(length(x) - tlag * (dimension-1)),
    n.reference = as.integer(n.reference),
    n.neighbor  = NA,
    metric      = metric,
    tlag        = as.integer(tlag),
    olag        = as.integer(olag),
    series.name = series.name,
    series      = asVector(x),
    xlab        = "ln(density)",
    ylab        = "< ln(neighborhood radius) >",
    invariant   = "information dimension")
}

###
# embedSeries
###

"embedSeries" <- function(x, dimension=NULL, tlag=NULL, series.name=NULL)
{
  if (is(x,"embedSeries"))
    return(x)

  # obtain series name
  if (is.null(series.name))
    series.name <- deparseText(substitute(x))

  if (is.matrix(x)){

    z <- x
    n.embed   <- nrow(x)
    dimension <- ncol(x)
    tlag      <- NA
  }
  else{

    if (is.null(dimension))
      dimension <- 2
    if (is.null(tlag))
      tlag <- timeLag(x, method="acfzero")

    checkScalarType(dimension,"integer")
    checkScalarType(tlag,"integer")
    checkScalarType(series.name,"character")

    x <- as.numeric(x)
    storage.mode(x) <- "double"

    N <- length(x)

    z <- itCall( "RS_fractal_embed",
      as.numeric(x), as.integer(dimension), as.integer(tlag))
      #COPY=rep(FALSE,3), #CLASSES=c("matrix", "integer", "integer"),
      #PACKAGE="ifultools")

    dimension <- ncol(z)
    n.embed   <- nrow(z)
    dim.names <- character(dimension)

    for (i in seq(length=dimension)){

      dim.names[i] <- switch(as.integer(i > 1) + 1,
			     paste(series.name,"[t]", sep=""),
			     paste(series.name,"[t + ", (i-1) * tlag[1], "]", sep=""))
    }

    dimnames(z) <- list(paste("t=",seq(n.embed),sep=""), dim.names)
  }

  oldClass(z) <- "embedSeries"

  attr(z, "emb.dim")     <- seq(dimension)
  attr(z, "n.embed")     <- n.embed
  attr(z, "tlag")        <- tlag
  attr(z, "series.name") <- series.name

  z
}

###
# FNN
###

"FNN" <- function(x, dimension=5, tlag=NULL, rtol=10, atol=2, olag=1)
{
  # obtain series name
  if (is(x,"signalSeries"))
    data.name <- x@title
  else
    data.name <- deparseText(substitute(x))

  # coerce to a vector
  x <- asVector(x)

  # check input arguments
  if (is.null(tlag))
    tlag <- timeLag(x, method="acfdecor")
  if (olag <= 0)
    stop("For FNN, the orbital lag must be a positive integer")

  # check input argument types and lengths
  checkScalarType(dimension,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(rtol,"numeric")
  checkScalarType(atol,"numeric")
  checkScalarType(olag,"integer")

  z <- itCall( "RS_fractal_dimension_false_nearest_neighbors",
    asVector(x), as.integer(dimension), as.integer(tlag), as.integer(olag), as.numeric(rtol), as.numeric(atol))
    #COPY=rep(FALSE,6),
    #CLASSES=c("matrix", rep("integer", 3), rep("numeric", 2)),
    #PACKAGE="ifultools")

  dimnames(z) <- list(c("rtol","atol","combined"), paste("E=",seq(dimension),sep=""))

  # assign class
  oldClass(z) <- "FNN"

  attr(z, "atol") <- atol
  attr(z, "rtol") <- rtol
  attr(z, "olag") <- olag
  attr(z, "tlag") <- tlag
  attr(z, "data.name") <- data.name
  attr(z, "n.sample")  <- length(x)
  attr(z, "dimension") <- dimension

  z
}

###
# FNN
###

"FNS" <- function(x, dimension=5, tlag=NULL, atol=1, image.tol=1, olag=1)
{
  # obtain series name
  if (class(x) == "signalSeries")
    data.name <- x@title
  else
    data.name <- deparseText(substitute(x))

  # coerce to a vector
  x <- asVector(x)

  # check input arguments
  if (is.null(tlag))
    tlag <- timeLag(x, method="acfdecor")
  if (olag <= 0)
    stop("For FNS, the orbital lag must be a positive integer")

  # check input argument types and lengths
  checkScalarType(dimension,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(atol,"numeric")
  checkScalarType(olag,"integer")
  checkScalarType(image.tol,"integer")
  if (image.tol	< 0)
    stop("image.tol must be non-negative")

  z <- itCall("RS_fractal_dimension_false_nearest_strands",
    as.numeric(x), as.integer(dimension), as.integer(tlag), as.integer(olag), as.integer(image.tol), as.numeric(atol))
    #COPY=rep(FALSE,6), #CLASSES = c("matrix", rep("integer", 4), "numeric"),
    #PACKAGE="ifultools")

  dimnames(z) <- list(paste("E=",seq(dimension),sep=""), "FNS %")

  # assign class
#  oldClass(z) <- "FNN"

#  attr(z, "atol") <- atol
#  attr(z, "rtol") <- rtol
#  attr(z, "olag") <- olag
#  attr(z, "tlag") <- tlag
#  attr(z, "data.name") <- data.name
#  attr(z, "n.sample")  <- length(x)
#  attr(z, "dimension") <- dimension

  z
}

###
# lyapunov
###

"lyapunov" <- function(x, tlag=NULL, dimension=5, local.dimension=3,
   reference=NULL, n.reference=NULL, olag=2,
   sampling.interval=NULL, polynomial.order=3, metric=Inf, scale=NULL)
{
  # obtain series name
  series.name <- deparseText(substitute(x))

  if (!is.numeric(x) && !is(x,"signalSeries"))
    stop ("Input data must be a numeric vector or an object of class \"signalSeries\"")

  # define defaults for missing inputs
  if (is.null(sampling.interval))
    sampling.interval <- deltat(x)
  if (is.null(tlag))
    tlag <- as.integer(timeLag(x, method="acfdecor"))
  if (is.null(n.reference))
    n.reference <- min(as.integer(round(length(x) / 20)), 100)

  # check input argument types and lengths
  checkScalarType(tlag,"integer")
  checkScalarType(dimension,"integer")
  checkScalarType(local.dimension,"integer")
  checkScalarType(olag,"integer")
  checkScalarType(sampling.interval,"numeric")
  checkScalarType(polynomial.order,"integer")
  checkScalarType(metric,"numeric")

  # check input arguments
  if (dimension < 1)
    stop("Embedding dimension must be at least unity")
  if (local.dimension < 1)
    stop("Local dimension must be at least unity")
  if (local.dimension > dimension)
    stop("Local dimension cannot exceed embedding dimension")
  if (polynomial.order < 1)
    stop("Polynomial order must be at least unity")
  if (n.reference < 10)
    stop("Number of reference points msut be at least 10")
  if (tlag < 1)
    stop("Time lag must be at least unity")
  if (sampling.interval <= 0.0)
    stop("Sampling interval must be a positive numeric value")

  # initialize variables
  n.embed   <- length(x) - (dimension - 1) * tlag
  scale.max <- n.embed - 2 - n.reference

  if (is.null(scale))
    scale <- as.integer(2^(seq(min(floor(logb(scale.max,2)) - 2 , 10)) - 1))

  scale <- as.integer(scale)

  if (any(scale < 1) || any(scale > scale.max))
    stop(paste("Scales must be on interval [1, ", scale.max, "]",sep=""))

  # ... reference (initial starting point in base index 1)
  if (is.null(reference))
    reference <- as.integer(
      seq(from=1, to=n.embed - max(scale) - n.reference - 2, length=5))

  reference <- as.integer(unique(reference))

  if (any(reference < 1))
	  stop("Initial reference points must be positive integers")
  if (any(reference > n.embed - max(scale) - n.reference + 1))
	  stop("Initial reference points are too large")
  if (olag < 1)
    stop("Orbital lag must be a positive integer")
  if (any(reference >= n.embed - max(scale)))
    stop("An initial reference point is set too large")

  # estimate the local Lyapunoiv spectrum about
  # each supplied reference
  z <- itCall( "RS_fractal_local_lyapunov_spectrum",
    asVector(x),
    as.integer(dimension),
    as.integer(tlag),
    as.integer(olag),
    as.numeric(sampling.interval),
    as.integer(local.dimension),
    as.integer(polynomial.order),
    as.integer(reference),
    as.integer(n.reference),
    mutilsDistanceMetric(metric),
    as.integer(scale))
    #COPY=rep(FALSE,11),
    #CLASSES=c("matrix", rep("integer", 3), "numeric", rep("integer", 5), "matrix"),
    #PACKAGE="ifultools")

  nms <- list(reference,scale)
  z <- lapply(z,
    function(x,nms){
	  x <- data.frame(x)
	  dimnames(x) <- nms
	  return(x)},
	nms=nms)

  oldClass(z) <- "lyapunov"
  attr(z, "scales")            <- names(z[[1]])
  attr(z, "n.exponent")        <- length(z)
  attr(z, "tlag")              <- tlag
  attr(z, "dimension")         <- dimension
  attr(z, "local.dimension")   <- local.dimension
  attr(z, "reference")         <- reference
  attr(z, "n.reference")       <- n.reference
  attr(z, "olag")              <- olag
  attr(z, "sampling.interval") <- sampling.interval
  attr(z, "polynomial.order")  <- polynomial.order
  attr(z, "metric")            <- metric
  attr(z, "n.sample")          <- length(x)
  attr(z, "data.name")         <- series.name

  return(z)
}

###
# poincareMap
###

"poincareMap" <- function(x, extrema="min", denoise=FALSE)
{
  checkScalarType(extrema,"character")

  # perform argument checks
  if (!(isVectorAtomic(x) && is.numeric(x)))
    stop("Input x must be a numeric vector")

  # map extrema
  extrema <- match.arg(extrema, c("min","max","all"))

  # calculate Poincare map via extrema
  z <- itCall( "RS_fractal_poincare_map",
    as.matrix(as.double(x)),
    as.integer(switch(extrema, min=0, max=1, all=2)),
    as.logical(denoise))
    #COPY=rep(FALSE,3),
    #CLASSES=c("matrix","integer","logical"),
    #PACKAGE="ifultools")

  z <- lapply(z, as.vector)
  names(z) <- c("location","amplitude")

  z
}

###
# spaceTime
###

"spaceTime" <- function(x, dimension=2, tlag=timeLag(x, method="acfdecor"),
  olag.max=as.integer(min(500,length(x)/20)), probability=0.1)
{
  # obtain series name
  if (is(x,"signalSeries"))
    series.name <- x@title
  else
    series.name <- deparseText(substitute(x))

  # ensure x is a numeric vector
  x <- asVector(x)
  storage.mode(x) <- "double"

  # perform argument checks
  checkScalarType(dimension,"integer")
  checkScalarType(tlag,"integer")
  checkScalarType(olag.max,"integer")
  checkScalarType(probability,"numeric")

  if (dimension < 1)
    stop ("Input dimension must be a positive integer")
  if (tlag < 1)
    stop("Input tlag must be a positive integer")
  if (olag.max < 1)
    stop("Input olag.max must be a positive integer")
  if ((probability <= 0) | (probability >= 1))
    stop("Input probability must be in the interval (0,1)")
  if ((numRows(x) - (dimension-1) * tlag) < 2)
    stop("Not enough embedding points")

  # initialize variables
  n.embed <- numRows(x) - (dimension-1) * tlag
  if (olag.max >= n.embed)
    olag.max <- n.embed - 1
  olags   <- seq(olag.max)

  # call wrapper function
  z <- t(itCall( "RS_fractal_space_time_separation_plot",
    as.matrix(x),
    as.integer(dimension),
    as.integer(tlag),
    as.matrix(as.integer(olags)),
    as.numeric(probability))[[2]])
    #COPY=rep(FALSE,5),
    #CLASSES=c("matrix","integer","integer","matrix","numeric"),
    #PACKAGE="ifultools")[[2]])

  oldClass(z) <- "spaceTime"

  attr(z, "dimension")   <- dimension
  attr(z, "n.embed")     <- n.embed
  attr(z, "tlag")        <- tlag
  attr(z, "probability") <- probability
  attr(z, "series.name") <- series.name
  attr(z, "olags")       <- olags

  z
}

###
# timeLag
###

"timeLag" <- function(x, method="acfzero", plot.data=FALSE)
{
  # local functions
  "firstACF" <- function(x, value=0){

    # calculate the autocorrelation function
    s      <- sapa::ACVS(x)
    data   <- s / s[1]
    lags   <- seq(0, length(data) - 1)
    icross <- which(data < value)[1]

    if (is.na(icross)){
      stop(paste("\n\nACF does not cross ", value,
        ". choose the time lag manually or select ",
		"another automated method.\n",sep=""))
    }
    else{

      # select closest point to given value
      if (abs(data[icross] - value) > abs(data[icross - 1] - value))
        icross <- icross - 1
    }

    return(list(closest=icross, data=data, lags=lags))
  }

  checkScalarType(method,"character")
  checkScalarType(plot.data,"logical")

  # obtain series name
  data.name <- deparseText(substitute(x))

  # coerce method string to lowercase
  method <- match.arg(lowerCase(method),
    c("acfzero","acfdecorrelation","acfnadir","mutual"))

  if (method == "acfzero"){

    obj  <- firstACF(x, value=0)
    z    <- obj$closest
    data <- obj$data
    lags <- obj$lags
  }
  else if (method == "acfdecorrelation"){

    obj  <- firstACF(x, value= 1 / exp(1))
    z    <- obj$closest
    data <- obj$data
    lags <- obj$lags
  }
  else if (method == "acfnadir"){

    # calculate the autocorrelation function
    s    <- sapa::ACVS(x)
    data <- s / s[1]
    lags <- seq(0, length(data) - 1)

    z <- which(as.logical(peaks(- data, strict=TRUE)))[1] - 1
  }
  else if (method == "mutual"){

    # calculate the mutual information for each lag from
    # 1:P where P is the 1.5 times the first nadir of the ACF
    high <- round(as.numeric(timeLag(x, method="acfnadir", plot.data=FALSE)) * 1.5)
    lags <- seq(1, high)

    data <- itCall( "RS_fractal_time_delayed_mutual_information", as.numeric(x), as.integer(lags))
      #COPY=rep(FALSE,2), #CLASSES=c("matrix", "matrix"),
      #PACKAGE="ifultools")
    z <- which(as.logical(peaks(-as.vector(data), strict=TRUE)))[1] - 1
  }

  # add attributes
  attr(z, "data")   <- data
  attr(z, "lags")   <- lags
  attr(z, "method") <- method

  # plot the results
  if (plot.data){

    if (charmatch("acf", method, nomatch=FALSE))
      ylab <- "Autocorrelation Function"
    else
      ylab <- "Time Delayed Mutual Information"

    cut <- seq(1, min(c(floor(1.5 * z), length(data))))

    plot(lags[cut], data[cut], type="h", xlab="lag", ylab=ylab)
    points(lags[cut], data[cut], pch=5)
    abline(h=0, lty=2)
    points(lags[z], data[z], pch=18, col=8, cex=2)
    title(paste("Best lag:", z))
  }

  z
}

