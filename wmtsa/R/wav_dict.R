################################################
## WMTSA wavelet dictionary functionality
##
##  Constructor Functions and methods:
##
##    wavDictionary
##
##      as.list.wavDictionary
##      print.wavDictionary
##
################################################

###
# wavDictionary (constructor)
###

"wavDictionary" <- function(wavelet, dual, decimate, n.sample,
  attr.x, n.levels, boundary, conv, filters, fast, is.complex)
{
  # check input arguments
  checkScalarType(wavelet,"character")
  checkScalarType(dual,"logical")
  checkScalarType(decimate,"logical")
  checkScalarType(n.sample,"integer")
  checkScalarType(n.levels,"integer")
  checkScalarType(boundary,"character")
  if (!is.list(filters))
    stop("filters must be a list")
  if (!all(is.element(c("scaling","wavelet"), names(filters))))
    stop("filters list must have named elements \"scaling\" and \"wavelet\"")
  if (!isVectorAtomic(filters$wavelet) || !isVectorAtomic(filters$scaling))
    stop("filters$wavelet and filters$scaling must be vectors ",
      "as defined by isVectorAtomic()")
  checkScalarType(fast,"logical")
  checkScalarType(is.complex,"logical")

  z <- list(wavelet=wavelet, dual=dual, decimate=decimate, n.sample=n.sample,
    attr.x=attr.x, n.levels=n.levels,
    boundary=boundary, conv=conv,
    analysis.filter=list(low=filters$scaling, high=filters$wavelet),
    synthesis.filter=list(low=filters$scaling, high=filters$wavelet),
    fast=fast, is.complex=is.complex)

  oldClass(z) <- "wavDictionary"

  return(z)
}

###
# as.list.wavDictionary
###

"as.list.wavDictionary" <- function(x,...)
{
	list(
	  "Wavelet"=x$wavelet,
	  "Length of series"=x$n.sample,
	  "Number of levels"=x$n.levels,
	  "Boundary correction rule"=x$boundary,
	  "Degree of polynomial"=ifelse1(x$boundary=="polynomial", attr(x$boundary, "pdeg"), NULL),
	  "Filtering technique"=ifelse1(x$conv,"convolution","correlation"))
}

###
# print.wavDictionary
###

"print.wavDictionary" <- function(x, ...)
{
  #define local functions

  "all.wavelet.names" <- function()
  {
    cnames <- c("c6","c12","c18","c24","c30")
    hnames <- c("haar","Haar")
    dnames <- c("d2", "d4","d6","d8","d10","d12","d14","d16","d18","d20")
    snames <- c("s4","s6", "s8","s10","s12","s14","s16","s18","s20")
    c(cnames, hnames, dnames, snames)
  }
  wn <- x$wavelet
  cat("Wavelet:", wn)
  if(!match(wn, all.wavelet.names(), nomatch=0)) print(x$analysis.filter, ...)
  cat("\nLength of series:", x$n.sample)
  cat("\nNumber of levels:", x$n.levels)
  boundary <- x$boundary
  cat("\nBoundary correction rule:", boundary)
  if(boundary=="polynomial")
    cat("\nDegree of polynomial:", attr(boundary, "pdeg"))

  # include filtering technique in printout
  conv <- x$conv
  cat("\nFiltering technique: ")
  if (conv) cat("convolution")
  else cat("correlation")
  cat("\n")

  invisible(x)
}



