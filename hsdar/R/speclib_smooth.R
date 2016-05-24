smoothSpeclib <- function(
                           x,
                           method="mean",
                           ...
                          )
{
  if (x@spectra@fromRaster)
    return(.blockwise(speclib_obj =  "x", pos = 1))
  
  predefinedmethod <- FALSE
  if (!is.speclib(x))
    stop("x must be of class 'Speclib'")
    
  if (!x@continuousdata)
    stop("Smoothing is only useful for continuous spectra")
  
  setmask <- FALSE  
  
  if (!is.null(attr(x, "setmask")))
  {
    if (attr(x, "setmask"))
    {
      setmask <- TRUE
      x <- interpolate.mask(x)
    }
  }
  
  res <- x
  
  if (method=="sgolay")
  {
    spectra(res) <- t(apply(spectra(x), 1, FUN = sgolayfilt, ...))
    usagehistory(res) <- paste("Smoothed with Savitzky-Golay smoothing filter")
    predefinedmethod <- TRUE
  }
  
  if (method=="lowess")
  {
    lowessFUN <- function(y, x, ...) lowess(x=x, y=y, ...)$y
    
    wavelength(res) <- lowess(x=wavelength(x), y=spectra(x)[1,], ...)$x
    spectra(res) <- t(apply(spectra(x), 1, FUN = lowessFUN, 
                            x=wavelength(x), ...))
    usagehistory(res) <- paste("Smoothed with lowess function")
    predefinedmethod <- TRUE
  } 
  if (method=="spline")
  {
    splineFUN <- function(y, x, ...) spline(x=x, y=y, ...)$y
    
    wavelength(res) <- spline(x=wavelength(x), y=spectra(x)[1,], ...)$x
    spectra(res) <- t(apply(spectra(x), 1, FUN = splineFUN, 
                            x=wavelength(x), ...))
    usagehistory(res) <- paste("Smoothed with spline function")
    predefinedmethod <- TRUE
  } 
  if (any(method==c("mean","mean_gliding")))
  {
    spectra  <- spectra(res)
    
    spectra(res) <- meanfilter(spectra, ...)
    usagehistory(res) <- paste("Smoothed with meanfilter")
    predefinedmethod <- TRUE
  }
  
  if (!predefinedmethod)
  {
    spectra(res) <- t(apply(spectra(x), 1, FUN = method, ...))
    usagehistory(res) <- paste("Smoothed with ", method," smoothing filter")
  }
  
  if (setmask) mask(res) <- attr(res, "dropped")
  return(res)            
}

meanfilter <- function(spectra, p=5)
{
  gliding  <- FALSE
  spectra  <- as.matrix(spectra)
  nwl      <- ncol(spectra)
  n        <- nrow(spectra)
  smoothed <- spectra*0
  
  storage.mode(nwl)      <- "integer"
  storage.mode(n)        <- "integer"
  storage.mode(p)        <- "integer"
  storage.mode(spectra)  <- "double"
  storage.mode(smoothed) <- "double"  

  if (!gliding)
  {
    external <- .Fortran("meanfilter",
                        nwl=nwl,
                        n=n,
                        p=p,
                        y=spectra,
                        smoothed=smoothed,
                        PACKAGE="hsdar"
                        )
  } else {
    external <- .Fortran("gliding_meanfilter",
                        nwl=nwl,
                        n=n,
                        p=p,
                        y=spectra,
                        smoothed=smoothed,
                        PACKAGE="hsdar"
                        )
  }
  
  external$smoothed <- as.data.frame(external$smoothed)
  
  return(external$smoothed)
}