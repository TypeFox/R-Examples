setMethod('writeStart', signature(x = 'HyperSpecRaster', filename = "character"),
           function(x, filename, ...) 
 {
   mc <- match.call(definition = sys.function(-1), 
                    call = sys.call(-1), expand.dots = TRUE)
   args <- as.list(mc)[-1]
   if (any(names(args) == "nl"))
   {
     x <- as(x, "RasterBrick")
     x <- brick(x, nl = eval(parse(text = args$nl)))
   }
   return(callNextMethod(x, filename, ...))
 } 
)

setMethod('writeValues', signature(x='HyperSpecRaster', v='Speclib'),
           function(x, v, start) 
 {
   if (nlayers(x) != nbands(v))
     stop(paste("Number of bands in x (", 
                nlayers(x), ") and v (",
                nbands(v),") differ", sep =""))
   spec <- spectra(v)   
   return(writeValues(x, spec, start))
 } 
)

setMethod('writeValues', signature(x='RasterBrick', v='Speclib'),
           function(x, v, start) 
 {
   if (nlayers(x) != nbands(v))
     stop(paste("Number of bands in x (", 
                nlayers(x), ") and v (",
                nbands(v),") differ", sep =""))
   spec <- spectra(v)
   return(writeValues(x, spec, start))
 } 
)

setMethod('writeValues', signature(x='RasterLayer', v='Speclib'),
           function(x, v, start) 
 {
   if (nbands(v) != 1)
     stop(paste("Number of bands in v must equal 1 if x is of class RasterLayer", sep =""))
   spec <- spectra(v)   
   return(writeValues(x, spec, start))
 } 
)

setMethod('show', 'HyperSpecRaster',
           function(object) 
 {
   callNextMethod(object)
   cat('wavelength:\n')
   print(object@wavelength)
 } 
)

setMethod('getValuesBlock', signature(x = 'HyperSpecRaster'),
           function(x, ...) 
 {
   v <- callNextMethod(x, ...)
   return(speclib(v, x@wavelength, fwhm = if (length(x@fwhm) > 0) x@fwhm else NULL, 
                  attributes = if (nrow(x@attributes) > 0) x@attributes else NULL)) 
 } 
)

setMethod('HyperSpecRaster', signature(x = 'character', wavelength = "numeric"),
           function(x, wavelength, fwhm = NULL, attributes = NULL, ...) 
 {
   res <- brick(x, ...)
   if (nlayers(res) != length(wavelength))
     stop("Length of wavelength do not equal number of bands in file")
   res <- as(res, 'HyperSpecRaster')
   res@wavelength <- wavelength
   if (!is.null(fwhm))
     res@fwhm <- fwhm
   if (!is.null(attributes))
     res@attributes <- attributes
   return(res)
 } 
)

setMethod('HyperSpecRaster', signature(x = 'RasterBrick', wavelength = "numeric"),
           function(x, wavelength, fwhm = NULL, attributes = NULL) 
 {
   if (nlayers(x) != length(wavelength))
     stop("Length of wavelength do not equal number of bands in file")
   res <- as(x, 'HyperSpecRaster')
   res@wavelength <- wavelength
   if (!is.null(fwhm))
     res@fwhm <- fwhm
   if (!is.null(attributes))
     res@attributes <- attributes
   return(res)
 } 
)

setMethod('HyperSpecRaster', signature(x = 'RasterLayer', wavelength = "numeric"),
           function(x, wavelength, fwhm = NULL, attributes = NULL) 
 {
   res <- brick(x, nl = length(wavelength))
   res <- as(res, 'HyperSpecRaster')
   res@wavelength <- wavelength
   if (!is.null(fwhm))
     res@fwhm <- fwhm
   if (!is.null(attributes))
     res@attributes <- attributes
   return(res)
 } 
)


setMethod('brick', signature(x = 'Speclib'),
           function(x, nrow, ncol, xmn, xmx, ymn, ymx, crs) 
 {
   if (missing(nrow))
   {
     if (!.is.rastermeta(x))
     {
       stop("nrow missing")
     } else {
       nrow <- x@rastermeta[[1]][1]
     }
   }
   if (missing(ncol))
   {
     if (!.is.rastermeta(x))
     {
       stop("ncol missing")
     } else {
       ncol <- x@rastermeta[[1]][2]
     }
   }
   if (missing(xmn))
   {
     if (!.is.rastermeta(x))
     {
       xmn <- 0 
     } else {
       xmn <- x@rastermeta[[2]]@xmin
     }
   }
   if (missing(xmx))
   {
     if (!.is.rastermeta(x))
     {
       xmx <- 1 
     } else {
       xmx <- x@rastermeta[[2]]@xmax
     }
   }
   if (missing(ymn))
   {
     if (!.is.rastermeta(x))
     {
       ymn <- 0 
     } else {
       ymn <- x@rastermeta[[2]]@ymin
     }
   }
   if (missing(ymx))
   {
     if (!.is.rastermeta(x))
     {
       ymx <- 1
     } else {
       ymx <- x@rastermeta[[2]]@ymax
     }
   }
   if (missing(crs))
   {
     if (!.is.rastermeta(x))
     {
       crs <- NA
     } else {
       crs <- x@rastermeta[[3]]
     }
   }
   
   
   arr <- array(data = spectra(x), dim = c(ncol, nrow, nbands(x)))
   res <- brick(arr, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs, transpose = TRUE)

   return(res)
 } 
)

setMethod('HyperSpecRaster', signature(x = 'Speclib'),
            function(x, nrow, ncol, xmn, xmx, ymn, ymx, crs) 
 {
   if (missing(nrow))
   {
     if (!.is.rastermeta(x))
     {
       stop("nrow missing")
     } else {
       nrow <- x@rastermeta[[1]][1]
     }
   }
   if (missing(ncol))
   {
     if (!.is.rastermeta(x))
     {
       stop("ncol missing")
     } else {
       ncol <- x@rastermeta[[1]][2]
     }
   }
   if (missing(xmn))
   {
     if (!.is.rastermeta(x))
     {
       xmn <- 0 
     } else {
       xmn <- x@rastermeta[[2]]@xmin
     }
   }
   if (missing(xmx))
   {
     if (!.is.rastermeta(x))
     {
       xmx <- 1 
     } else {
       xmx <- x@rastermeta[[2]]@xmax
     }
   }
   if (missing(ymn))
   {
     if (!.is.rastermeta(x))
     {
       ymn <- 0 
     } else {
       ymn <- x@rastermeta[[2]]@ymin
     }
   }
   if (missing(ymx))
   {
     if (!.is.rastermeta(x))
     {
       ymx <- 1
     } else {
       ymx <- x@rastermeta[[2]]@ymax
     }
   }
   if (missing(crs))
   {
     if (!.is.rastermeta(x))
     {
       crs <- NA
     } else {
       crs <- x@rastermeta[[3]]
     }
   }
   
   
   arr <- array(data = spectra(x), dim = c(ncol, nrow, nbands(x)))
   res <- brick(arr, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs, transpose = TRUE)
   wl <- wavelength(x)
   if (is.data.frame(wl))
     wl <- rowMeans(wl[,c(1,2)])
      
   return(HyperSpecRaster(res, wavelength = wl))
 }
)


setReplaceMethod("wavelength", signature(object = "HyperSpecRaster", value = "numeric"), 
                 function(object, value)
{
  object@wavelength <- value
  return(object)
}
)

setMethod("wavelength", signature(object = "HyperSpecRaster"), 
          function(object)
  return(object@wavelength)
)