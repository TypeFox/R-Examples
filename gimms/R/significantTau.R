if ( !isGeneric("significantTau") ) {
  setGeneric("significantTau", function(x, ...)
    standardGeneric("significantTau"))
}
#' Compute (pre-whitened) Kendall's tau
#'
#' @description
#' Apply the Mann-Kendall trend test (Mann, 1945) to a series of observations
#' and return Kendall's tau (Kendall, 1938) based on a predefined significance
#' level. In contrast to other readily available implementations, it is left to
#' the user to decide whether or not to apply pre-whitening as described in the
#' \strong{zyp} package vignette (Bronaugh and Werner, 2013).
#'
#' @param x A 'numeric' vector.
#' @param p Significance level to be tested.
#' @param prewhitening Logical. If \code{TRUE}, pre-whitening is applied prior
#' to the Mann-Kendall trend test.
#' @param df Logical. If \code{TRUE}, a 'data.frame' holding the value of
#' Kendall's tau and the referring significance level.
#' @param filename Character. Optional output filename. Needs to include an
#' appropriate file format extension, see \code{\link{writeFormats}}.
#' @param ... Further arguments passed on to \code{\link{zyp.trend.vector}}.
#'
#' @return If \code{df = FALSE} (default) and \code{p} was not exceeded, a
#' single 'numeric'; if \code{df = FALSE} and \code{p} was exceeded, a 'logical'
#' (\code{NA}); else a 'data.frame' with Kendall's tau and the corresponding
#' significance level.
#'
#' @author
#' Florian Detsch
#'
#' @seealso
#' \code{\link{MannKendall}}, \code{\link{zyp.trend.vector}}.
#'
#' @references
#' Kendall, M.G. (1938). A new measure of rank correlation. \emph{Biometrika}
#' 30(1/2), 81-93, doi: 10.2307/2332226. Available online at
#' \url{http://www.jstor.org/stable/2332226} (accessed 2015-11-06).
#'
#' Mann, H.B. (1945). Nonparametric tests against trend. \emph{Econometrica}
#' 13(3), 245-259, doi: 10.2307/1907187. Available online at
#' \url{http://www.jstor.org/stable/1907187} (accessed 2015-11-06).
#'
#' Zhang, X., Vincent, L.A., Hogg, W.D. and A. Niitsoo (2000). Temperature and
#' Precipitation Trends in Canada during the 20th Century. \emph{Atmosphere-Ocean}
#' 38(3), 395-429, doi: 10.1080/07055900.2000.9649654. Available online at
#' \url{http://www.tandfonline.com/doi/abs/10.1080/07055900.2000.9649654}
#' (accessed 2015-11-06).
#'
#' Yue, S., Pilon, P., Phinney, B. and G. Cavadias (2002). The influence of
#' autocorrelation on the ability to detect trend in hydrological series.
#' \emph{Hydrological Processes} 16, 1807-1829, doi: 10.1002/hyp.1095.
#' Available online at \url{http://onlinelibrary.wiley.com/doi/10.1002/hyp.1095/abstract}
#' (accessed 2015-11-06).
#'
#' @examples
#' ## Example taken from ?Kendall::MannKendall
#' library(Kendall)
#' data(PrecipGL)
#' plot(PrecipGL)
#'
#' ## Mann-Kendall trend test without pre-whitening
#' x <- as.numeric(PrecipGL)
#' significantTau(x, p = 0.001, prewhitening = FALSE, df = TRUE)
#'
#' ## Mann-Kendall trend test with pre-whitening
#' significantTau(x, p = 0.001, prewhitening = TRUE, df = TRUE)
#'
#' #############################################################################
#' ### use case: significant (p < 0.001) mann-kendall trends (2009-13) #########
#' #############################################################################
#'
#' \dontrun{
#' ## download files from 2009-2013
#' gimms_files <- downloadGimms(x = as.Date("2009-01-01"),
#'                              dsn = paste0(getwd(), "/data"))
#'
#' ## convert binary files to 'Raster*' format
#' gimms_rasters <- rasterizeGimms(gimms_files, filename = gimms_files,
#'                                 format = "GTiff", overwrite = TRUE)
#'
#' ## crop iran
#' library(rworldmap)
#' data(countriesLow)
#' gimms_iran <- crop(gimms_rasters,
#'                      subset(countriesLow, ADMIN == "Iran"))
#'
#' ## remove seasonal signal via remote::deseason
#' library(remote)
#' gimms_deseason <- deseason(gimms_iran, cycle.window = 24, use.cpp = TRUE)
#'
#' ## identify long-term monotonous trends from mann-kendall trend test; all values
#' ## of kendall's tau with p >= 0.001 are set to NA; note that pre-whitening is
#' ## applied prior to the actual trend test
#' gimms_trends <- significantTau(gimms_deseason, p = 0.001,
#'                                prewhitening = TRUE)
#'
#' ## create figure
#' library(RColorBrewer)
#' library(latticeExtra)
#' cols <- colorRampPalette(brewer.pal(11, "BrBG"))
#' spplot(gimms_trends, col.regions = cols(100), scales = list(draw = TRUE),
#'        at = seq(-.525, .525, .025),
#'        main = expression("Kendall's" ~ tau ~ "(2009-2013)")) +
#'   layer(sp.polygons(countriesLow))
#' }
#'
#' @export significantTau
#' @name significantTau

################################################################################
### function using 'numeric' input #############################################
#' @aliases significantTau,numeric-method
#' @rdname significantTau
setMethod("significantTau",
          signature(x = "numeric"),
          function(x, p = 0.001, prewhitening = TRUE, df = FALSE, ...) {

            # if only one unique value exists in 'x', return NA
            if (length(unique(x)) == 1)
              return(NA)

            # with prewhitening
            if (prewhitening) {

              # try to compute pre-whitened mann-kendall trend test
              try(mk <- zyp::zyp.trend.vector(x, ...), silent = TRUE)

              # if previous computation fails, return NA
              if (!exists("mk")) {
                sig <- tau <- NA
                # else return kendall's tau and referring p value
              } else {

                id_sig <- grep("sig", names(mk))
                sig <- mk[id_sig]

                id_tau <- grep("tau", names(mk))
                tau <- mk[id_tau]
              }

              # without prewhitening
            } else {
              mk <- Kendall::MannKendall(x)

              sig <- mk$sl
              tau <- mk$tau
            }

            # return data.frame
            if (df) {
              return(data.frame(tau = tau, p = sig))

              # reject value of tau if p >= 0.001
            } else {

              if (is.logical(sig) | is.logical(p)) {
                return(NA)
              } else {
                if (sig >= p) {
                  return(NA)
                  # keep value of tau if p < 0.001
                } else {
                  return(tau)
                }
              }
            }
          }
)


################################################################################
### function using 'RasterStack' or 'RasterBrick' input ########################
#' @aliases significantTau,RasterStackBrick-method
#' @rdname significantTau
setMethod("significantTau",
          signature(x = "RasterStackBrick"),
          function(x, p = 0.001, prewhitening = TRUE, filename = "", ...) {

            rst_mk <- raster::overlay(x, fun = function(y, ...) {
              significantTau(y, p = p,
                             prewhitening = prewhitening, df = FALSE, ...)
            })

            ## write to file
            if (nchar(filename) > 0)
              rst_mk <- raster::writeRaster(rst_mk, filename = filename,
                                            overwrite = TRUE)

            ## return
            return(rst_mk)
          }
)
