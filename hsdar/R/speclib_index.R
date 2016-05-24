setMethod("[", "Speclib",
          function(x, i, j, ...)
{
  dots <- list(...)
  upduh <- !any(names(dots) == "usagehistory")
  if (missing(i)) 
  {
    tmp <- spectra(x, j = j)
    if (nspectra(x) == 1)
    { 
      x@spectra@fromRaster <- FALSE      
      spectra(x) <- matrix(tmp, ncol = length(tmp))
    } else {
      x@spectra@fromRaster <- FALSE
      spectra(x) <- tmp
    }
    wavelength(x) <- wavelength(x)[j]

    if (upduh)
      usagehistory(x) <- "Subsetting speclib (spectral dimension)"
    return(x)
  } 
  if (missing(j))
  {
    tmp <- spectra(x, i = i) 
    if (class(tmp) == "numeric")
      tmp <- matrix(tmp, ncol = if (nbands(x) > 1) length(tmp) else 1)
    x@spectra@fromRaster <- FALSE
    spectra(x) <- tmp
    idSpeclib(x) <- as.character(idSpeclib(x)[i])
    at_x <- attribute(x)[i,]
    if (! class(at_x) %in% c("matrix", "data.frame"))
    {
      at_x <- data.frame(x = at_x)
      names(at_x) <- names(attribute(x))
    }
    attribute(x) <- at_x   

    if (upduh)
      usagehistory(x) <- "Subsetting speclib (sample dimension)"
    return(x)
  }
  tmp <- spectra(x, i = i, j = j)
  if (class(tmp) == "numeric")
  {
    ncols <- sum(rep.int(1, nbands(x))[j])
    nrows <- sum(rep.int(1, nspectra(x))[i])
    tmp <- matrix(tmp, ncol = ncols, nrow = nrows)
  }
  x@spectra@fromRaster <- FALSE
  spectra(x) <- tmp
  wavelength(x) <- wavelength(x)[j]
  idSpeclib(x) <- as.character(idSpeclib(x)[i])

  at_x <- attribute(x)[i,]
  if (! class(at_x) %in% c("matrix", "data.frame"))
  {
    at_x <- data.frame(x = at_x)
    names(at_x) <- names(attribute(x))
  }
  attribute(x) <- at_x   

  if (upduh)
    usagehistory(x) <- "Subsetting speclib (spectral and sample dimensions)"
  return(x)
})