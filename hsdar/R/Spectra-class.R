setMethod("ncol", signature(x = ".Spectra"), 
          function(x)
  if (x@fromRaster)
  {
    return(x@spectra_ra@data@nlayers)
  } else {
    return(ncol(x@spectra_ma))
  }
)

setMethod("nrow", signature(x = ".Spectra"), 
          function(x)
  if (x@fromRaster)
  {
    return(x@spectra_ra@nrows * x@spectra_ra@ncols)
  } else {
    return(nrow(x@spectra_ma))
  }
)