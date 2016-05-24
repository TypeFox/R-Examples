setMethod("merge", signature(x = "Speclib", y = "Speclib"),
          function(x, y, ...)
{
  if (dim(x)[2] != dim(y)[2])
    stop("Dimensions of Speclibs do not fit")

  wl <- wavelength(x)
  if (any(wl!=wavelength(y)))
    stop("Wavelengths differ")

  if (nrow(attribute(y)) == dim(y)[1])
  {
    if (nrow(attribute(x)) == dim(x)[1])
    {
      attribute(x) <- rbind(attribute(x),attribute(y))
    } else {
      warning("x does not have proper attributes definition. Attributes information will be lost")
    }
  } else {
    warning("y does not have proper attributes definition. Attributes information will be lost")
  }
  
  spectra(x) <- as.matrix(rbind(spectra(x),spectra(y))) 
  
  return(x)
}
)