##' Center and scale hyperSpec object
##'
##' \code{link[base]{scale}}s the spectra matrix. \code{scale (x, scale = FALSE)} centers the data.
##'
##' Package \code{scale} provides a fast alternative for \code{base::\link[base]{scale}}
##' 
##' @name scale,hyperSpec-method
##' @rdname scale
##' @aliases scale scale-methods scale,hyperSpec-method
##' @docType methods
##' @param x the \code{hyperSpec} object
##' @param center if \code{TRUE}, the data is centered to \code{colMeans (x)}, \code{FALSE}
##' suppresses centering. Alternatively, an object that can be converted to numeric of length
##' \code{nwl (x)} by \code{\link[base]{as.matrix}} (e.g. hyperSpec object containing 1 spectrum) can
##' specify the center spectrum.
##' @param scale if \code{TRUE}, the data is scaled to have unit variance at each wavelength,
##' \code{FALSE} suppresses scaling. Alternatively, an object that can be converted to numeric of
##' length \code{nwl (x)} by \code{\link[base]{as.matrix}} (e.g. hyperSpec object containing 1 spectrum)
##' can specify the center spectrum.
##' @return the centered & scaled \code{hyperSpec} object
##' @author C. Beleites
##' @seealso \code{\link[base]{scale}}
##'
##' package scale.
##' @keywords methods
##' @export
##' @examples
##'
##' ## mean center & variance scale
##' tmp <- scale (chondro)
##' plot (tmp, "spcmeansd")
##' plot (sample (tmp, 5), add = TRUE, col = 2)
##'
##' ## mean center only
##' tmp <- scale (chondro, scale = FALSE)
##' plot (tmp, "spcmeansd")
##' plot (sample (tmp, 5), add = TRUE, col = 2)
##'
##' ## custom center
##' tmp <- sweep (chondro, 1, mean, `/`)
##' plot (tmp, "spcmeansd")
##' tmp <- scale (tmp, center = quantile (tmp, .05), scale = FALSE)
##' 
setMethod ("scale", signature = signature (x = "hyperSpec"),
           function (x, center = TRUE, scale = TRUE){
  validObject (x)

  if (! is.logical (center)) center <- as.matrix (center)
  if (! is.logical (scale))  scale  <- as.matrix (scale)
  
  x@data$spc <- scale (x@data$spc, center, scale)

  x
})
