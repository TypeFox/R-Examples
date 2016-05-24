setClass(".Spectra",
         representation(
           fromRaster = "logical",
           spectra_ma = "matrix",
           spectra_ra = 'RasterBrick'
         ),
         prototype(
           fromRaster = FALSE,
           spectra_ma = matrix(),
           spectra_ra = new("RasterBrick")
         )
)

setClass("Speclib",
         representation(
           spectra = ".Spectra", 
           wavelength = "numeric",
           attributes = "data.frame",
           fwhm = "numeric",
           continuousdata = "logical",
           wlunit = "character",
           xlabel = "character",
           ylabel = "character",
           ID = "character",
           wavelength.is.range = "logical",
           transformation = "character",
           usagehistory = "character",
           rastermeta = "list"
         ),
         prototype(
           spectra = new(".Spectra"),                   
           wavelength = numeric(),
           attributes = data.frame(),
           fwhm = 1,
           continuousdata = TRUE,
           wlunit = "nm",
           xlabel = "Wavelength",
           ylabel = "Reflectance",
           ID = character(),
           wavelength.is.range = FALSE,
           transformation = "NONE",
           usagehistory = "",
           rastermeta = list()
         ),
         validity = function(object)
         {
           c1 <- TRUE
           if (nrow(object@attributes) > 0)
           {
             c1 <- nrow(object@attributes) == nrow(object@spectra)
             if (!c1)
             {
               stop("Invalid attribute data.frame for spectra")
             }
           }
           c2 <- ncol(object@spectra) == length(object@wavelength)
           if (!c2)
           {
             stop("Invalid wavelength vector for spectra")
           }
           return(c1 & c2)
         }
)

setClass('HyperSpecRaster',
         contains = 'RasterBrick',
         representation(
           wavelength = 'numeric',
           fwhm       = 'numeric',
           attributes = 'data.frame'
         ),
         prototype (
           wavelength = numeric(),
           fwhm       = numeric(),
           attributes = data.frame()
                   )
        )

setClass("DistMat3D",
         representation(
           values = "numeric",
           ncol = "numeric",
           nlyr = "numeric"
         ),
         prototype(
           values = numeric(),
           ncol = 0,
           nlyr = 0
         ),
         validity = function(object)
         {
           if (length(object@values) != object@nlyr * (sum(1:object@ncol)-object@ncol))
             stop("Length of values do not fit dimensions of matrix")
         }
)

setClass("Nri",
         representation(
           nri = "DistMat3D",
           fwhm = "numeric",
           wavelength = "numeric",
           dimnames = "list",
           multivariate = "list",
           attributes = "data.frame",
           usagehistory = "character"
         ),
         prototype(
           nri = new("DistMat3D", values = numeric(), nlyr = 0, ncol = 0),
           fwhm = 0,
           wavelength = 0,
           dimnames = list(),
           multivariate = list(),                   
           attributes = data.frame(),
           usagehistory = ""
         ),
         validity = function(object)
         {
           if (length(object@wavelength) != object@nri@ncol)
             stop("Length of wavelength do not fit dimensions of nri")
         }
)

setClass('Clman',
         contains = 'Speclib',
         representation(
           cp         = 'matrix',
           hull       = 'matrix'
         ),
         prototype (
           cp         = matrix(),
           hull       = matrix()
                   ),
         validity = function(object)
         {
           if (ncol(object@cp) != length(object@wavelength))
             stop("Number of bands in continuum points and length of wavelength differ")
           if (ncol(object@spectra) != length(object@wavelength))
             stop("Number of bands in spectra and length of wavelength differ")
           if (nrow(object@spectra) != nrow(object@cp))
             stop("Number of samples in spectra and continuum points differ") 
           return(TRUE)
         }
        )

setClass('Specfeat',
         contains = 'Speclib',
         representation(
           features      = 'list',
           featureLimits = 'list'
         ),
         prototype (
           features      = list(),
           featureLimits = list()
                   )#,
#          validity = function(object)
#          {
#            if (ncol(object@cp) != length(object@wavelength))
#              stop("Number of bands in continuum points and length of wavelength differ")
#            if (ncol(object@spectra) != length(object@wavelength))
#              stop("Number of bands in spectra and length of wavelength differ")
#            if (nrow(object@spectra) != nrow(object@cp))
#              stop("Number of samples in spectra and continuum points differ") 
#            return(TRUE)
#          }
        )

setClassUnion(".CaretHyperspectral", c("Speclib", "Nri", "Specfeat"))

if (!isGeneric("speclib")) {
  setGeneric("speclib", function(spectra, wavelength, ...)
  standardGeneric("speclib"))
}

if (!isGeneric("spectra")) {
  setGeneric("spectra", function(object, ...)
  standardGeneric("spectra"))
}
if (!isGeneric("spectra<-")) {
  setGeneric("spectra<-",function(object, value)
  standardGeneric("spectra<-"))
}

if (!isGeneric("mask")) {
  setGeneric("mask", function(object, ...)
  standardGeneric("mask"))
}
if (!isGeneric("mask<-")) {
  setGeneric("mask<-",function(object, value)
  standardGeneric("mask<-"))
}

if (!isGeneric("attribute")) {
  setGeneric("attribute", function(object)
  standardGeneric("attribute"))
}
if (!isGeneric("attribute<-")) {
  setGeneric("attribute<-", function(object, value)
  standardGeneric("attribute<-"))
}

if (!isGeneric("wavelength")) {
  setGeneric("wavelength", function(object, ...)
  standardGeneric("wavelength"))
}
if (!isGeneric("wavelength<-")) {
  setGeneric("wavelength<-",function(object, value)
  standardGeneric("wavelength<-"))
}

if (!isGeneric("distMat3D")) {
  setGeneric("distMat3D",function(x, ...)
  standardGeneric("distMat3D"))
}

if (!isGeneric('HyperSpecRaster')) 
{
  setGeneric('HyperSpecRaster', function(x, wavelength, ...)
  standardGeneric('HyperSpecRaster')) 
}


if (!isGeneric("ncol")) {
  setGeneric("ncol", function(object, ...)
  standardGeneric("ncol"))
}

if (!isGeneric("nrow")) {
  setGeneric("nrow", function(object, ...)
  standardGeneric("nrow"))
}

if (!isGeneric("as.data.frame")) {
  setGeneric("as.data.frame")
}