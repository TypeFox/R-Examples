setMethod("speclib", signature(spectra = "matrix", wavelength = "numeric"), 
          function(spectra, wavelength, ...)
  return(createspeclib(spectra, wavelength, ...))
)

setMethod("speclib", signature(spectra = "SpatialGridDataFrame", wavelength = "numeric"), 
          function(spectra, wavelength, ...)
  {
    spectra <- t(as.matrix(.getImgMatrix_SpatialGridDataFrame(spectra)))
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "numeric", wavelength = "numeric"), 
          function(spectra, wavelength, ...)
  {
    spectra <- matrix(spectra, ncol = length(wavelength))
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "matrix", wavelength = "data.frame"), 
          function(spectra, wavelength, ...)
  return(createspeclib(spectra, wavelength, ...))
)

setMethod("speclib", signature(spectra = "SpatialGridDataFrame", wavelength = "data.frame"), 
          function(spectra, wavelength, ...)
  {
    spectra <- t(as.matrix(.getImgMatrix_SpatialGridDataFrame(spectra)))
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "numeric", wavelength = "data.frame"), 
          function(spectra, wavelength, ...)
  {
    spectra <- matrix(spectra, ncol = nrow(wavelength))
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "matrix", wavelength = "matrix"), 
          function(spectra, wavelength, ...)
  return(createspeclib(spectra, wavelength, ...))
)

setMethod("speclib", signature(spectra = "SpatialGridDataFrame", wavelength = "matrix"), 
          function(spectra, wavelength, ...)
  {
    spectra <- t(as.matrix(.getImgMatrix_SpatialGridDataFrame(spectra)))
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "numeric", wavelength = "matrix"), 
          function(spectra, wavelength, ...)
  {
    spectra <- matrix(spectra, ncol = nrow(wavelength))
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "HyperSpecRaster"), 
          function(spectra, ...)
  {
    ref_system <- as.character(crs(spectra))
    
    v <- getValues(spectra)
    
    res <- speclib(v, spectra@wavelength, fwhm = if (length(spectra@fwhm) > 0) spectra@fwhm else NULL, 
                   attributes = if (nrow(spectra@attributes) > 0) spectra@attributes else NULL, 
                   rastermeta = rastermeta(spectra))
    return(res)                                    
  }
)

setMethod("speclib", signature(spectra = "RasterBrick", wavelength = "numeric"), 
          function(spectra, wavelength, ...)
  {
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "RasterBrick", wavelength = "data.frame"), 
          function(spectra, wavelength, ...)
  {
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("speclib", signature(spectra = "RasterBrick", wavelength = "matrix"), 
          function(spectra, wavelength, ...)
  {
    return(createspeclib(spectra, wavelength, ...))
  }
)

setMethod("$", signature(x = "Speclib"), 
          function(x, name)
  {
    slot(x, name)
  }
)

setMethod("as.data.frame", signature(x = "Speclib"), 
          function(x, row.names = NULL, optional = FALSE, includeAttr = FALSE, ...)
  {
    x_dat <- as.data.frame(spectra(x), row.names = row.names, optional = optional, ...)
    names(x_dat) <- bandnames(x)
    if (!includeAttr)
      return(x_dat)
      
    if (ncol(attribute(x)) > 0)
      return(cbind(x_dat, attribute(x)))
  }
)


createspeclib <- function (spectra,
                           wavelength,
                           fwhm = NULL,
                           attributes = NULL,
                           usagehistory = NULL,
                           transformation = NULL,
                           continuousdata = "auto",
                           wlunit = "nm",
                           xlabel = "Wavelength",
                           ylabel = "Reflectance",
                           rastermeta = NULL
                          )
{
  if (class(spectra) %in% c("RasterBrick", "HyperSpecRaster"))
  {
    fromRaster <- TRUE
  } else {
    fromRaster <- FALSE
    dim_spectra <- c(nrow(spectra), ncol(spectra))
  }
  wavelength.is.range <- FALSE
  if (class(wavelength)=="data.frame" || class(wavelength)=="matrix")
  {
    wavelength <- as.data.frame(wavelength)
    if (ncol(wavelength)==1)
    {
      wavelength <- as.vector(wavelength)
      if (!is.null(fwhm))
        if (length(fwhm)!=length(wavelength))
          stop("Length of fwhm and wavelength differ")
      if (!fromRaster)
      {
        if (length(wavelength)==dim_spectra[1])
        {
          if (length(wavelength)==dim_spectra[2])
          {
            warning("Could not determine orientation of spectra data. \n  
                    Make sure that columns are wavelength and rows samples")
          } else {
            spectra <- t(spectra)
          }
        }
      }
    } else {
      wavelength.is.range <- TRUE
      if (!fromRaster)
      {
        if (nrow(wavelength)==dim_spectra[1])
        {
          if (nrow(wavelength)==dim_spectra[2])
          {
            warning("Could not determine orientation of spectra data. \n  
                    Make sure that columns are wavelength and rows samples")
          } else {
            spectra <- t(spectra)
          }
        }
      }
    }
  } else {
    if (!is.null(fwhm))
    {
      wavelength.is.range <- TRUE
      if (length(fwhm)!=length(wavelength))
        stop("Length of fwhm and wavelength differ")
    }
    if (!fromRaster)
    {
      if (length(wavelength)==dim_spectra[1])
      {
        if (length(wavelength)==dim_spectra[2])
        {
          warning("Could not determine orientation of spectra data. \n  
                  Make sure that columns are bands and rows different samples")
        } else {
          spectra <- t(spectra)
        }
      }
    }
  }
  if (!fromRaster)
  {
    names <- NULL
    rn <- row.names(spectra)
    if (!is.null(rn))
    {
      rn <- as.factor(rn)
      if (nlevels(rn)!=length(rn))
      {
        warning("  some row.names duplicated: --> Spectra does not have IDs")
        rn <- as.factor(1:length(rn))
      }
      rn <- as.character(rn)
    } else {
      rn <- character()
    }
    spectra <- as.matrix(spectra)
    cn <- colnames(spectra)
    rownames(spectra) <- NULL
    colnames(spectra) <- NULL
  } else {
    cn <- c(1:spectra@data@nlayers)
    rn <- character()    
  }
  
  
  if (!wavelength.is.range)
  {
    range <- wavelength[-1] - wavelength[-1*length(wavelength)]
    range <- c(as.numeric(range),range[length(range)])
    if (sd(range)==0)
      range <- mean(range)
    fwhm <- range
  } else {
    if (!is.null(fwhm))
    {
      if (is.data.frame(wavelength))
        wavelength <- rowMeans(wavelength)
    }
  }
    
  if (is.null(attributes)) 
    attributes <- data.frame()
    
  if (is.null(usagehistory))   
    usagehistory <- character()
  
  if (is.null(transformation))   
    transformation <- character()
    
  if (is.null(rastermeta))   
    rastermeta <- list()
  
  result <- new("Speclib", 
                spectra = spectra, 
                wavelength = wavelength,
                fwhm = fwhm,
                wavelength.is.range = wavelength.is.range,
                continuousdata = continuousdata,
                attributes = attributes,
                transformation = transformation,
                usagehistory = usagehistory,
                wlunit = wlunit,
                xlabel = xlabel,
                ylabel = ylabel,
                rastermeta = rastermeta                
               )
  idSpeclib(result) <- rn
  bandnames(result) <- cn
  if (validObject(result))
  { 
    return(result)
  }
}

setMethod("initialize", signature(.Object = "Speclib"),
          function(.Object, ...)
{  
  dots <- list(...)
  if (any(names(dots) == "continuousdata"))
  {
    if (dots$continuousdata != "auto")
    {
      if (mode(dots$continuousdata) != "logical")
        stop("continuousdata must be 'auto', TRUE or FALSE")
      continuousdata <- dots$continuousdata
    } else {
      continuousdata <- max(dots$wavelength[-1*length(dots$wavelength)]-dots$wavelength[-1]) <= 20
    }
  } else {
    continuousdata <- TRUE
  }
  
  if (any(names(dots) == "wavelength"))
  {
    wavelength <- dots$wavelength
  } else {
    wavelength <- numeric()
#     stop("Wavelength information required")
  }

  if (any(names(dots) == "spectra"))
  {
    spectra <- dots$spectra
    fromRaster <- class(spectra) %in% c("RasterBrick", "HyperSpecRaster")
    spectra <- new(".Spectra",
                   fromRaster = fromRaster,
                   spectra_ma = if (fromRaster) matrix() else spectra,
                   spectra_ra = if (fromRaster) spectra else new("RasterBrick"))
  } else {
    spectra <- new(".Spectra")
#     stop("Spectra required")
  }

  if (any(names(dots) == "attributes"))
  {
    attributes <- dots$attributes
  } else {
    attributes <- data.frame()
  }

  if (any(names(dots) == "fwhm"))
  {
    fwhm <- dots$fwhm
  } else {
    fwhm <- dots$wavelength[-1] - dots$wavelength[-1*length(dots$wavelength)]
    fwhm <- c(fwhm, fwhm[length(fwhm)])
  }

  if (any(names(dots) == "wavelength.is.range"))
  {
    wavelength.is.range <- dots$wavelength.is.range
  } else {
    wavelength.is.range <- FALSE
  }

  if (any(names(dots) == "transformation"))
  {
    transformation <- dots$transformation
  } else {
    transformation <- "NONE"
  }
  
  if (any(names(dots) == "usagehistory"))
  {
    usagehistory <- dots$usagehistory
  } else {
    usagehistory <- ""
  }

  if (any(names(dots) == "wlunit"))
  {
    wlunit <- dots$wlunit
  } else {
    wlunit <- "nm"
  }

  if (any(names(dots) == "xlabel"))
  {
    xlabel <- dots$xlabel
  } else {
    xlabel <- "Wavelength"
  }

  if (any(names(dots) == "ylabel"))
  {
    ylabel <- dots$ylabel
  } else {
    ylabel <- "Reflectance"
  }
  
  if (any(names(dots) == "rastermeta"))
  {
    rastermeta <- dots$rastermeta
  } else {
    rastermeta <- list()
  }
  
  object <- .Object
  
  object@spectra             <- spectra
  object@wavelength          <- wavelength
  object@fwhm                <- fwhm
  object@continuousdata      <- continuousdata
  object@wavelength.is.range <- wavelength.is.range
  object@transformation      <- transformation
  object@attributes          <- attributes
  object@usagehistory        <- usagehistory
  object@wlunit              <- wlunit
  object@xlabel              <- xlabel
  object@ylabel              <- ylabel
  object@rastermeta          <- rastermeta
  return(object)
}
)

.getImgMatrix_SpatialGridDataFrame <- function(sp_dat)
  return(as.matrix(sp_dat@data))
  
  

