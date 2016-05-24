# dist.default <- dist

# dist <- function(x,...) UseMethod("dist")

dist.speclib <- function(
                         x,
                         method="sam",
                         ...
                         )
{
  if (class(x)!="Speclib") 
    stop("x must be of class 'Speclib'")
  
  if (method=="sam")
  {
    distance <- sam_distance(x)
    distance <- as.dist(distance)
  } else {    
    if (attr(x, "setmask"))
      x <- interpolate.mask(x)
      
    spec <- spectra(x)
    
    distance <- dist(spec, method = method, ...)
  }
  return(distance)
}


sam <- function(
                x,
                ref
               )
{
  if (x@spectra@fromRaster)
    return(.blockwise(speclib_obj =  "x", pos = 1))
  
  if (class(x)!="Speclib") 
    stop("x must be of class 'Speclib'")
  if (class(ref)!="Speclib")
    stop("ref must be of class 'Speclib'")
    
  spec <- spectra(x)
  wlx <- x@wavelength
  
  
  specref <- spectra(ref)
  wlref <- ref@wavelength

  if (length(wlref) != length(wlx))
  {
    stop("Wavelength between speclibs differ")
  }
  
  spec    <- as.matrix(spec)
  specref <- as.matrix(specref)

  if (max(spec, na.rm = TRUE)>1)
  {
    spec <- spec/100
    specref <- specref/100
  }
  
  if (max(spec, na.rm = TRUE)>1)
    stop("Spectra in x must be in range [0,1]")
  if (max(specref, na.rm = TRUE)>1)
    stop("Spectra in ref must be in range [0,1]")
    
  nspec   <- nrow(spec)
  nref    <- nrow(specref)
  nbands  <- ncol(spec)
  specang <- array(0, dim = c(nspec,nref))
  

  storage.mode(nspec)    <- "integer"
  storage.mode(nref)     <- "integer"
  storage.mode(nbands)   <- "integer"
  storage.mode(spec)     <- "double"
  storage.mode(specref)  <- "double"
  storage.mode(specang)  <- "double"
  
  distance <- .Fortran("sam",
                       nspec=nspec,
                       nref=nref,
                       nbands=nbands,
                       spec=spec,
                       specref=specref,
                       specang=specang,
                       PACKAGE="hsdar"
                       )$specang
 
  distance <- as.matrix(distance)
  colnames(distance) <- rownames(specref)
  rownames(distance) <- rownames(spec)
  return(distance)                             
}

sam_distance <- function (x)
{
  if (class(x)!="Speclib") 
    stop("x must be of class 'Speclib'")
  
    
  spec <- spectra(x)
  
  if (attr(x, "setmask"))
    x <- interpolate.mask(x)
  
  spec    <- as.matrix(spec)  
  nspec   <- nrow(spec)
  nbands  <- ncol(spec)
  specang <- array(0, dim = c(nspec,nspec))
  if (max(spec)>1)
    spec <- spec/100
  
  storage.mode(nspec)    <- "integer"
  storage.mode(nbands)   <- "integer"
  storage.mode(spec)     <- "double"
  storage.mode(specang)  <- "double"
    
  distance <- .Fortran("sam",
                       nspec=nspec,
                       nref=nspec,
                       nbands=nbands,
                       spec=spec,
                       specref=spec,
                       specang=specang,
                       PACKAGE="hsdar"
                       )$specang
 
  distance <- as.matrix(distance)
  colnames(distance) <- rownames(spec)
  rownames(distance) <- rownames(spec)
  diag(distance) <- 0
  return(distance)     
}
