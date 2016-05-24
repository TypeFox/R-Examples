################################################
## WMTSA utility functionality
##
##   oceansdf
##   wavMaxLevel
##   wavSortCrystals
##   wavTitle
##
################################################

###
# oceansdf
###

"oceansdf" <- function(f)
{
  ## spectral density function for ocean shear data

  ## develop frequency ranges
  ff <- abs(f)
  f.low <- f[ff <= (5/128)]
  f.mid <- f[ff > (5/128) & ff <= (5/8)]
  f.high <- f[ff > (5/8) & ff <= 5]

  ## develop spectral density function
  S.low <- rep(32.768, length(f.low))
  S.mid <- 32.768 * abs((128 * f.mid)/5)^(-3.5)
  S.high <- (abs(8/5 * f.high)^(-1.7))/500
  S <- c(S.low, S.mid, S.high)
  if(length(S) < length(f))
    warning("Some frequencies out of range")
  S
}

###
# wavMaxLevel
###

"wavMaxLevel" <- function(n.taps=8, n.sample=1024, xform="modwt")
{
  xform <- match.arg(lowerCase(xform), c("modwt","modwpt","dwt","dwpt"))
  xform <- ifelse1(xform == "dwpt", "dwt", xform == "modwpt", "modwt", xform)
  type  <- mutilsTransformType(xform)

  n.interior <- itCall("RS_wavelets_transform_coefficient_boundaries",
     as.integer(20), as.integer(n.taps), as.integer(n.sample), as.integer(type))[[3]]
    #
     #COPY=rep(FALSE,4),
     #CLASSES=rep("integer",4),
     #PACKAGE="ifultools")[[3]]

   # here we require that the number of interior wavelet coefficients
   # be greater than unity (rather than greater than zero). we do this
   # because it leads to problems in the C code for wavelet variance in
   # the case that the n.coefficients is unity.

   as.integer(max(which(as.logical(n.interior > 1))))
}

###
# wavSortCrystals
###

"wavSortCrystals" <- function(x, reverse=FALSE)
{
  # check for proper class
  if (!is.character(x))
    stop("Input must be of class \"character\"")

  # develop local sort function
  localsort <- function(x)
    x[order(as.numeric(substring(x,2)))]

  # divide the string into wavelet and scaling crystals
  wave.crystals <- x[grep("d",x)]
  scal.crystals <- x[grep("s",x)]

  # now sort each individually and recombine
  sorted.crystals <- c(localsort(wave.crystals), localsort(scal.crystals))

  # reverse if requested
  if (!reverse)
    return(sorted.crystals)
  else
    return(rev(sorted.crystals))
}

###
# wavTitle
###

"wavTitle" <- function(x, default="x")
{
  # define local functions
  "wavStripString" <- function(x)
  {
    # strip the result by keeping the inner most
    # name within parentheses. For example,
    # "a(b(c(d, other_args)))" would return "d".
    # this assumes that the first argument
    # always contains the input time series argument.

    # turn string into a character vector
    # each element containing one character
    nch <- nchar(x)
    str <- substring(x, 1:nch, 1:nch)

    # now find the location of the innermost
    # "(" character
    str.start <- which(as.logical(charmatch(str,"(")))

    if (all(!is.missing(str.start))){

      str.start <- max(str.start)

      # remove first portion of string up to
      # the innermost "("
      str <- str[(min(str.start+1,nch): nch)]

      # redefine string length
      nch <- length(str)

      # now test for the end of the argument
      # which will be either where the next "," or
      # ")" occurs, whichever comes first

      right.paren <- min(which(as.logical(charmatch(str,")"))),1)
      right.comma <- min(which(as.logical(charmatch(str,","))),1)

      # find which occurs first, the comma or the right
      # parenthesis
      str.end <- min(c(right.paren, right.comma), na.rm=TRUE)

      # return the stripped substring
      if (!is.na(str.end))
        y <- paste(str[(1: max(1, str.end - 1))], collapse="")
      else
        y <- paste(str, collpase="")
    }
    else
      y <- x

    y
  }

  # obtain class of input data
  wave.class <- class(x)

  # based on the class, attempt to extract the
  # name of the data
  if (is.element(wave.class, c("wavBoundary",
	"wavVarTest","wavTransform"))){
    series.name <- x$series@title
    alternative <- x$dictionary$series.name
  }
  else if (is.element(wave.class, "signalSeries")){
    series.name <- x@title
    alternative <- default
  }
  else if (is.element(wave.class, "wavFDP")){
    series.name <- x$series@title
    alternative <- default
  }
  else
    series.name <- deparseText(substitute(x))

  # convert the title
  if (length(series.name) < 1)
    series.name <- alternative

  if (is.null(series.name))
    series.name <- default

  # strip the series name
  #series.name <- wavStripString(series.name)

  series.name
}

