# if (!isGeneric("get_reflectance")) {
#   setGeneric("get_reflectance", function(spectra, ...)
#   standardGeneric("get_reflectance"))
# }

get_reflectance <- function(spectra, wavelength, position, weighted = FALSE, ...)
{
  if (wavelength[1]<=position & wavelength[length(wavelength)]>=position)
  {
    if (weighted)
    {
      if (any(wavelength==position))
      {
        return(get_reflectance(spectra, wavelength, position, weighted = FALSE))
      } else {
        temp <- abs(wavelength-position)
        ord <- order(temp)
        return((spectra[,ord[1]]*1/temp[ord[1]]+spectra[,ord[2]]*1/temp[ord[2]])/
               (1/temp[ord[1]]+1/temp[ord[2]]))
      }
    } else {
      temp <- abs(wavelength-position)
      return(spectra[,which(temp==min(temp))])
    }
  } else {
    return(rep.int(NA,nrow(spectra)))
  }
}

setMethod("get_reflectance", signature(spectra = "Speclib"), 
          function(spectra, position, ...)
{
  wavelength <- if (is.data.frame(spectra@wavelength)) rowMeans(spectra@wavelength) else spectra@wavelength
  spectra <- spectra(spectra)
  return(get_reflectance(spectra, wavelength, position, ...))  
}
)