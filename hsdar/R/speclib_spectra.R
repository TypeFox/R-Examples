setMethod("spectra", signature(object = "Speclib"), 
          function(object, i, j, ...)
{
  if (all(c(missing(i), missing(j))))
  {
    if (object@spectra@fromRaster)
    {
      if (object@spectra@spectra_ra@data@inmemory)
      {
        return(.spectra(object, ...)) 
      } else {
        return(object@spectra[])  
      }
    } else {
      return(.spectra(object, ...))            
    }
  } else {    
    if (missing(j))
      j <- c(1:nbands(object))
    if (object@spectra@fromRaster)
    {
      if (object@spectra@spectra_ra@data@inmemory)
      {
        if (missing(i))
          i <- c(1:nspectra(object))
        return(.spectra(object, ...)[i,j]) 
      } else {
        if (missing(i))
        {
          return(object@spectra[,j])
        } else {
          i <- c(1:nspectra(object))
          return(object@spectra[i,j])
        }
      }
    } else {
      if (missing(i))
        i <- c(1:nspectra(object))
      return(.spectra(object, ...)[i,j])            
    }
  }  
}
)

setReplaceMethod("spectra", signature(object = "Speclib", value = "matrix"), 
                 function(object, value)
{
  if (object@spectra@fromRaster)
  {
    object@spectra@spectra_ra <- setValues(object@spectra@spectra_ra, value)
  } else {
    object@spectra@spectra_ma <- value
  }
  return(object)
}
)

setReplaceMethod("spectra", signature(object = "Speclib", value = "data.frame"),
                 function(object, value)
{
  if (object@spectra@fromRaster)
  {
    object@spectra@spectra_ra <- setValues(object@spectra@spectra_ra, as.matrix(value))
  } else {
    object@spectra@spectra_ma <- as.matrix(value)
  }
  return(object)
}
)

setReplaceMethod("spectra", signature(object = "Speclib", value = "numeric"),
                 function(object, value)
{
  if (object@spectra@fromRaster)
  {
    object@spectra@spectra_ra <- setValues(object@spectra@spectra_ra,
                                           matrix(value, ncol = nbands(object)))
  } else {
    object@spectra@spectra_ma <- as.matrix(value)
  }
  return(object)
}
)


setReplaceMethod("spectra", signature(object = "Speclib", value = "RasterBrick"),
                 function(object, value)
{
  object@spectra@fromRaster <- TRUE
  object@spectra@spectra_ra <- value  
  return(object)
}
)


.spectra <- function(object, return_names = FALSE)
{
  if (object@spectra@fromRaster)
  {
    spec <- getValues(object@spectra@spectra_ra)
  } else {
    spec <- object@spectra@spectra_ma
  }

  if (return_names) 
  {    
    if (!is.null(bandnames(object)))
    {
      if (length(bandnames(object)) == ncol(spec))
      {
        colnames(spec) <- bandnames(object)
      } else {
        warning("Length of bandnames and number of bands in spectra differ. Drop bandnames")
        colnames(spec) <- paste("B_", wavelength(object), sep = "")
      }
    } else {
      colnames(spec) <- paste("B_", wavelength(object), sep = "")
    }
    
    if (length(idSpeclib(object)) > 0)
    {
      if (length(idSpeclib(object)) == nrow(spec))
      {
        rownames(spec) <- idSpeclib(object)
      } else {
        warning("Length of idSpeclib and number of spectra differ. Drop ID")
        rownames(spec) <- paste("ID_", c(1:nspectra(object)), sep = "")
      }
    } else {
      rownames(spec) <- paste("ID_", c(1:nspectra(object)), sep = "")
    }
  }
  return(spec)
}

setMethod("[", signature(x = ".Spectra"), 
          function(x, i, j, ...)
{
  if (!x@fromRaster)
    return(callNextMethod(x@spectra_ma, i, j, ...)) 
    
  if (missing(j)) 
  {
    if (missing(i)) 
      return(getValues(x@spectra_ra))
    j <- c(1:nlayers(x@spectra_ra))
  }
  if (is.logical(j))
    j <- c(1:length(j))[j]
    
  if (missing(i)) 
  {
    row <- 1
    nrows <- nrow(x@spectra_ra)
    col <- 1
    ncols <- ncol(x@spectra_ra)
    res <- getValuesBlock(x@spectra_ra, row = row, nrows = nrows, col = col, 
                          ncols = ncols, lyrs = j)
    
  } else {
    if (is.logical(i))
      i <- c(1:length(i))[i]
    
    cr <- rowColFromCell(x@spectra_ra, i)
    rows <- as.numeric(levels(as.factor(cr[,1])))
    ncols <- ncol(x@spectra_ra)
    res <- apply(matrix(rows, ncol = 1), 1, function(row, x, cols, lyr, ncols)
    {
      cols <- cols[cols[,1] == row,2]
      if (length(cols) > 1)
      {
        res <- getValuesBlock(x, row=row, nrows=1, col=1, ncols=ncols, lyr=lyr)
        res <- t(res[cols,])
      } else {
        res <- getValuesBlock(x, row=row, nrows=1, col=cols, ncols=1, lyr=lyr)
      }
      return(res)
    }, x@spectra_ra, cr, j, ncols)
    if (length(rows) == 1)
    {
      res <- matrix(res, ncol = length(j), byrow = TRUE)
    } else {
      if (is.matrix(res))
      {
        res <- t(res)
      } else {   
        res <- matrix(unlist(res), ncol = length(j), byrow = TRUE)
      }
    }
#     res <- t(apply(cr, 1, function(cr, x, lyr)
#     {
#       return(getValuesBlock(x, row=cr[1], nrows=1, col=cr[2], ncols=1, lyr=lyr))
#     }, x@spectra_ra, j))
  }
  return(res)     
})

setMethod("print", signature(x = ".Spectra"), 
          function(x)
{
  show(x)
}
)

setMethod("show", signature(object = ".Spectra"), 
          function(object)
{
  if (object@fromRaster)
  {
    cat("Spectra stored in RasterBrick object\n")
    print(object@spectra_ra)
    cat("Use 'spectra(x)[]' to read all data. Be careful if data is large!\n")
  } else {
    cat("Spectra stored in memory as matrix\n")
    print(object@spectra_ma)
  }
})