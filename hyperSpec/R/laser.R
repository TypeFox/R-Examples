##' Laser Emission
##' A time series of an unstable laser emission.
##' 
##' see the Vignette
##' 
##' @name laser
##' @docType data
##' @format The data set consists of 84 laser emission spectra measured during
##'   95 min.  Each spectrum has 36 data points in the range of 404.5 nm to
##'   405.8 nm.
##' @author C. Beleites
##' @keywords datasets
##' @examples
##' 
##' laser
##' 
##' cols <- c ("black", "blue", "darkgreen", "red")
##' wl <- c (405.0, 405.1, 405.3, 405.4)
##' plotspc (laser, axis.args=list (x = list (at = seq (404.5, 405.8, .1))))
##' for (i in seq_along (wl))
##'    abline (v = wl[i], col = cols[i], lwd = 2, lty = 2)
##' 
##' plotc (laser [,, wl], spc ~ t, groups = .wavelength, type = "b",
##'        col = cols)
##' 
##' \dontrun{vignette ("laser", package="hyperSpec")}
##' 
NULL
