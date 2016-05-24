################################################
## IFULTOOLS package MUTILS functions
##
##   mutilsDistanceMetric
##   mutilsFilterType
##   mutilsFilterTypeContinuous
##   mutilsSDF
##   mutilsTransformPeakType
##   mutilsTransformType
##   mutilsWaveletVarianceSDF
##   mutilsWSShrinkageFunction
##   mutilsWSThresholdFunction
##
###################################################

###
# mutilsDistanceMetric
##

"mutilsDistanceMetric" <- function(metric){

  checkScalarType(metric,"numeric")
  if (metric >= 1 && metric <= 2)
    metric <- metric - 1
  else if (metric == Inf)
    metric <- 2
  else
    stop("Acceptable values for metric are: 1, 2, and Inf")

  as.integer(metric)
}

###
# mutilsFilterType
###

"mutilsFilterType" <- function(wavelet="s8")
{
  # parse wavelet filter specs
  if (lowerCase(wavelet) == "haar"){
    filter.length <- 2
    filter.type   <- 7
  }
  else {
    # check for the correct filter types
    supported <- c("d","s","l","c")
    prefix    <- substring(wavelet,1,1)
    if (!is.element(prefix, supported))
      stop("The specified filter is currently unsupported")

    # map filter to index and obtain length
    filter.type   <- match(prefix, supported) - 1
    filter.length <- as.integer(substring(wavelet,2))
  }

  return(list(type=filter.type, length=filter.length))
}

###
# mutilsFilterTypeContinuous
###

"mutilsFilterTypeContinuous" <- function(x)
{
	x <- match.arg(lowerCase(x),c("haar","gauss1","gauss2",
	  "gaussian1","gaussian2","sombrero","mexican hat","morlet"))

  filter <- switch(x,
    haar=7,
    gauss1=4,
    gaussian1=4,
    gaussian2=5,
    gauss2=5,
    sombrero=5,
    "mexican hat"=5,
    morlet=6,
    NULL)

  as.integer(filter)
}

###
# mutilsSDF
###

"mutilsSDF" <- function(sdf=NULL, sdfargs=NULL, n.freq=1024, sampling.interval=1)
{
  # calculate sdf function over appropriate range of frequencies
  # such that f=[0, 1/P , 2/P, 3/P, ..., (M-1)/P] where P=2*(M-1)
  # (M = n.freq as defined above)
  if (!is.null(sdf) && is(sdf,"function")){

    checkScalarType(n.freq,"integer")
    checkScalarType(sampling.interval,"numeric")
    checkRange(n.freq, range.=c(2,Inf))
    checkRange(sampling.interval, range.=c(0,Inf), inclusion=c(FALSE,TRUE))
    if (!is.null(sdfargs) && !is.list(sdfargs))
      stop("sdfargs must either be NULL or a list of named additional inputs for the supplied SDF function")

    # check that "f" is the first argument of the input SDF function
    if (names(as.list(args(sdf)))[[1]] != "f")
      stop("Input SDF function must have as its first argument the variable \"f\"")

    f  <- seq(from=0, to=1/2, length=n.freq) / sampling.interval
    tmpargs <- c(list(f=f), sdfargs)
    Sx <- do.call("sdf", tmpargs)

    attr(Sx, "frequency") <- f
  }
  else{

    # this is a flag for the MUTILS C code indicating that no SDF is supplied.
    # it is convenient to specify a negative value since we assume all
    # input SDFs are positive valued
    Sx <- -1.0
  }

  Sx
}

###
# mutilsTransformPeakType
###

"mutilsTransformPeakType" <- function(x)
{
	checkScalarType(x,"character")
	supported.xforms <- c("extrema","maxima","minima")
	x <- match.arg(x, supported.xforms)
	match(x, supported.xforms) - 1
}

###
# mutilsTransformType
###

"mutilsTransformType" <- function(x)
{
	checkScalarType(x,"character")
	supported.xforms <- c("modwt","modwpt","dwt","dwpt")
	x <- match.arg(x, supported.xforms)
	match(x, supported.xforms) - 1
}

###
# mutilsWSShrinkageFunction
###

"mutilsWSShrinkageFunction" <- function(x)
{
	checkScalarType(x,"character")
	supported <- c("hard","soft","mid")
	x <- match.arg(x, supported)
	index <- as.integer(pmatch(x, supported) - 1)
	list(index=index, shrinkfun=x)
}

###
# mutilsWSThresholdFunction
###

"mutilsWSThresholdFunction" <- function(x)
{
	checkScalarType(x,"character")
	supported <- c("universal","minimax","adaptive")
	x <- match.arg(x, supported)
	index <- as.integer(pmatch(x, supported) - 1)
	list(index=index, threshfun=x)
}
