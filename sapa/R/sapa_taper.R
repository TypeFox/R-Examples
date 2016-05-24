################################################
## SAPA tapers and windows
##
##::::::::::::::::::::::::::::::::::::::::::::::
##
## Class: taper
## Constructor function: taper
## Methods:
##
##   as.matrix.taper
##   plot.taper
##   print.taper
##
################################################

###
# taper
###

"taper" <- function(type="rectangle", n.sample=100, n.taper=NULL, sigma=0.3, beta=4*pi*(n.sample-1)/n.sample,
  cutoff=floor(n.sample / 2), sidelobedB=80,
  roughness=n.sample/2, flatness=0.3, bandwidth=4, normalize=TRUE)
{
  # check input argument types and lengths
  checkScalarType(type,"character")
  checkScalarType(n.sample,"integer")
  checkScalarType(sigma,"numeric")
  checkScalarType(beta,"numeric")
  checkScalarType(cutoff,"numeric")
  checkScalarType(sidelobedB,"numeric")
  checkScalarType(roughness,"numeric")
  checkScalarType(flatness,"numeric")
  checkScalarType(bandwidth,"numeric")
  checkScalarType(normalize,"logical")

  # define local functions
  "tapersDPSS" <- function(N, n.tapers=2, NW=4)
  {
    N <- as.integer(N)
    if(!(N > 0))
      stop(paste("parameter N =", N, "should be positive integer"))
    n.tapers <- as.integer(n.tapers)
    if(!(n.tapers > 0))
      stop(paste("parameter n.tapers =", n.tapers, "should be positive integer"))
    if(n.tapers > N)
      stop(paste("parameter n.tapers =", n.tapers, "cannot be greater than\nparameter N =", N))
    if(n.tapers > 1) 
    {
        result <- .Call("R_sapa_dpss", as.integer(N), as.integer(n.taper), as.double(bandwidth), PACKAGE = "sapa")
        dim(result) <- c(N, n.tapers)
        result <- result[, n.tapers:1]
        for(k in 1:as.integer((n.tapers + 1)/2) * 2 - 1)
            if(!(sum(result[, k]) > 0)) result[, k] <-  - result[, k]
        tmp <- N - 1 - 2 * (0:(N - 1))
        for(k in 1:as.integer(n.tapers/2) * 2)
            if (!(sum(result[, k] * tmp) > 0))
        result[, k] <-  -result[, k]
    }
    else 
    {
        result <- tapersDPSS(N, n.taper=2, NW=NW)[,1]
        dim(result) <- c(N, 1)
    }
    result
  }

  # match taper type
  type <- lowerCase(type)

  # supported types: keep this order as it is commensurate
  # with the C code
  supported <- c("rectangle", "triangle","raised cosine", "hanning",
    "hamming", "blackman", "nuttall", "gaussian",
    "kaiser", "chebyshev", "born jordan", "sine",
    "parzen", "papoulis", "daniell", "dpss")

  taperstr <- match.arg(type, supported)

  # is it a multitaper?
  is.multitaper <- charmatch(taperstr, c("sine","dpss"), nomatch=FALSE)

  # define default for number fo tapers
  if (is.null(n.taper))
    n.taper <- as.integer(ifelse1(is.multitaper, 5, 1))
  checkScalarType(n.taper,"integer")

  # check multitapers
  if (!is.multitaper && n.taper > 1)
    stop("Multitapers only supported with sine and dpss tapers")


  # map additional taper parameter
  param <- as.double(0.0)

  # ... raised cosine
  if (taperstr == "raised cosine"){

    if (flatness < 0.0 || flatness > 1.0)
      stop("Fraction of flatness must be on [0,1] for raised cosine taper")

    param <- as.double(1.0 - flatness)
  }

  # ... Gaussian
  if (taperstr == "gaussian"){

    if (sigma <= 0.0)
      stop("Standard deviation of Gaussian must be positive.")

    param <- as.double(1.0 / sigma)
  }

  # ... Kaiser
  if (taperstr == "kaiser"){

    if (beta < 0.0)
      stop("Kaiser shape factor must be positive or equal to zero.")

    param <- as.double(beta)
  }

  # ... Chebyshev
  if (taperstr == "chebyshev"){

    if (sidelobedB <= 0.0)
      stop("Chebyshev sidelobedB bandwidth (in decibels) must be positive.")

    param <- as.double(sidelobedB)
  }

  # ... Parzen
  if (taperstr == "parzen"){

    cutoff <- as.integer(cutoff)

    if (cutoff <= 1)
      stop("Parzen cutoff must be greater than unity.")

    param <- as.double(cutoff)
  }

  # ... Papoulis
  if (taperstr == "papoulis"){

    cutoff <- as.integer(cutoff)

    if (cutoff <= 1)
      stop("Papoulis cutoff must be greater than unity.")

    param <- as.double(cutoff)
  }


  # ... Daniell
  if (taperstr == "daniell"){

    if (roughness < 0.0)
      stop("Daniell roughness must be positive.")

    param <- as.double(roughness)
  }

  if (n.taper < 1)
    stop("Number of tapers must be positive")
  if (n.sample < 1)
    stop("Number of samples must be positive")

  if (taperstr == "dpss"){

    z <- t(tapersDPSS(n.sample, n.taper=n.taper, NW=bandwidth))
  }
  else{

    z <- itCall("RS_signal_taper",
      match(taperstr,supported) - 1, n.taper, n.sample, param, normalize)
      #,
      #COPY=rep(FALSE,5),
      #CLASSES=c(rep("integer",3), "numeric", "logical"),
      #PACKAGE="ifultools")
  }

  if (n.taper == 1)
    z <- as.vector(z)

  oldClass(z) <- "taper"

  attr(z, "type")            <- taperstr
  attr(z, "n.sample")        <- n.sample
  attr(z, "n.taper")         <- n.taper
  attr(z, "gaussian.sigma")  <- sigma
  attr(z, "kaiser.beta")     <- beta
  attr(z, "normalize")       <- normalize
  attr(z, "chebyshev.sidelobedB") <- sidelobedB
  attr(z, "parzen.cutoff")   <- cutoff
  attr(z, "papoulis.cutoff") <- cutoff
  attr(z, "daniell.roughness") <- roughness
  attr(z, "dpss.bandwidth")  <- bandwidth
  attr(z, "cosine.flatness") <- flatness

  return(z)
}

###
# as.matrix.taper
###

"as.matrix.taper" <- function(x,...){
  xatt <- attributes(x)
  matrix(as.vector(x),nrow=xatt$n.taper,ncol=xatt$n.sample,...)
}

###
# plot.taper
###

"plot.taper" <- function(x, type="l", ylab=upperCase(attr(x,"type")), ...)
{
  xatt    <- attributes(x)
  n.taper <- xatt$n.taper

  checkScalarType(type,"character")
  checkScalarType(ylab,"character")

  if (n.taper > 1){

    data <- as.matrix(x)
    ylim <- range(data)
    lwd  <- seq(0.5,2,len=n.taper)

    for (j in seq(n.taper)){

      if (j == 1)
        plot(data[j,], ylim=ylim, type=type,lwd=lwd[j], ylab=ylab,...)
      else
        lines(data[j,], lwd=lwd[j],...)
    }
  }
  else{
    data <- as.vector(x)
    plot(data, type=type, ylab=ylab, ...)
  }

  invisible(NULL)
}

###
# print.taper
###

"print.taper" <- function(x, pre="", ...){

  xatt <- attributes(x)

  type <- switch(lowerCase(xatt$type),
    parzen="Lag Window",
    papoulis="Lag Window",
    "Taper")

  cat(type, ": ", xatt$type, "\n", sep="")
  cat(pre, "Number of points:", xatt$n.sample,"\n")
  cat(pre, "Number of tapers:", xatt$n.taper,"\n")
  cat(pre, "Normalized:", xatt$normalize,"\n")

  extra <- switch(lowerCase(xatt$type),
    gaussian=paste("Standard deviation:",xatt$gaussian.sigma),
    kaiser=paste("Shape parameter (beta):",xatt$kaiser.beta),
    chebyshev=paste("Sidelobe attenuation (dB):",xatt$chebyshev.sidelobedB),
    parzen=paste("Cutoff:",xatt$parzen.cutoff),
    papoulis=paste("Cutoff:",xatt$papoulis.cutoff),
    daniell=paste("Roughness factor:",xatt$daniell.roughness),
    dpss=paste("Bandwidth:",xatt$dpss.bandwidth),
    "raised.cosine"=paste("Flatness fraction:",xatt$cosine.flatness),
    "")

  if (nchar(extra) > 0)
    cat(pre, extra,"\n")

  invisible(x)
}
