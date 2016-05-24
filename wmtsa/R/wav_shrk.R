################################################
## WMTSA package wavelet shrinkage functionality
##
##  Functions:
##
##    wavShrink
##
################################################

###
# wavShrink
###

"wavShrink" <- function(x, wavelet="s8",
  n.level=ilogb(length(x), base=2),
  shrink.fun="hard", thresh.fun="universal", threshold=NULL,
  thresh.scale=1, xform="modwt", noise.variance=-1.0,
  reflect=TRUE)
{
  # initialize variables
  series.name <- deparse(substitute(x))

  # check input arguments
  if (is.element(class(x),c("named","signalSeries")))
    x <- as.vector(x)
  if (is.complex(x))
    stop("Time series must be real-valued")
  checkVectorType(x,"numeric")
  if (length(x) < 2)
    stop("Time series must contain more than one point")
  if (any(is.na(x)))
    stop("Time series contains NA values")
  checkScalarType(xform, "character")
  checkScalarType(thresh.fun, "character")
  checkScalarType(shrink.fun, "character")
  checkScalarType(wavelet, "character")
  checkScalarType(thresh.scale, "numeric")

  if (is.null(noise.variance))
    noise.variance <- -1.0 # flag to C-code to auto estimate noise variance
  checkScalarType(noise.variance, "numeric")

  xform <- match.arg(lowerCase(xform), c("dwt","modwt"))
  if (xform == "modwt"){
    thresh.fun <- "universal"
    decimated <- FALSE
  }
  else
    decimated <- TRUE

  if (n.level < 1)
    stop("Number of wavelet transform decomposition levels must be positive")
  if (n.level > ilogb(length(x), base=2))
    stop("Number of wavelet transform decomposition levels exceeds maximum")

  # obtain the wavelet and scaling filters. the wavelet argument
  # is checked here
  filters <- wavDaubechies(wavelet=wavelet, normalized=!decimated)[c("wavelet","scaling")]

  # initialize length parameters
  L <- length(filters$wavelet)
  N <- length(x)

  shrfun <- mutilsWSShrinkageFunction(shrink.fun)
  thrfun <- mutilsWSThresholdFunction(thresh.fun)

  if (thresh.scale <= 0)
    stop("Threshold scaling factor must be positive")

  if (!is.null(threshold)){

    # threshold values are explicitly set.
    # ensure that there are as many thresholds
    # as there are wavelet transform decomposition levels.
    # if a single threshold is given, replicate out to n.level.
    # and normalize appropriately (in the case of the MODWT this means
    # dividing the given threshold by 2^((j-1)/2) for level j
    checkVectorType(threshold, "numeric")
    if (length(threshold) != n.level){
      threshold <- ifelse1(xform == "modwt", threshold[1] / 2^(seq(0, n.level-1)/2),
        rep(threshold[1], n.level))
    }

    if (any(threshold <= 0))
      stop("All thresholds must be positive")
  }
  else{

  	# setting threshold to a single negative value
  	# prompts the C code to automatically
    # estimate the thresholds in the .Call
    threshold <- -1.0
  }

  if (reflect){

    # reflect the last Lj points of the series, where
    # Lj = (2^n.level - 1)(L - 1) + 1, and L=length(filters$wavelet)
    Lj <- (2^n.level - 1) * (L - 1) + 1
    ix <- seq(N-1, max(N-Lj,1))
    x  <- c(x, x[ix])
  }

  z <- as.vector(itCall("RS_wavelets_shrink", as.numeric(x), filters,
    as.double(matrix(threshold)),
    thrfun$index, as.numeric(thresh.scale), as.numeric(noise.variance),
    shrfun$index, as.integer(n.level), as.logical(decimated)))
    #
    #COPY=rep(FALSE,9),
    #CLASSES = c("matrix", "list", "matrix", "integer", rep("numeric",2),
	#    rep("integer",2), "logical"),
    #PACKAGE="ifultools"))

  if (reflect)
    z <- z[seq(N)]

  # assign attributes
  attr(z,"wavelet")        <- wavelet
  attr(z,"n.level")        <- n.level
  attr(z,"shrink.fun")     <- shrfun$shrinkfun
  attr(z,"thresh.fun")     <- thrfun$threshfun
  attr(z,"thresh.scale")   <- thresh.scale
  attr(z,"transform")      <- xform
  attr(z,"noise.variance") <- noise.variance
  attr(z,"reflect")        <- reflect

  z
}

