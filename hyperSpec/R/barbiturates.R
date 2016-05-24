##' Barbiturates Spectra from .spc example files
##' A time series of mass spectra in a list of hyperSpec objects.
##' 
##' 
##' @name barbiturates
##' @docType data
##' @format The data sets consists of 286 spectra. They are the result of importing the
##' BARBITUATES.SPC example data from Thermo Galactic's spc file format specification.
##' @author C. Beleites and Thermo Galactic
##' @references The raw data is available at \url{http://hyperspec.r-forge.r-project.org/blob/fileio.zip}
##' @keywords datasets
##' @examples
##' 
##' barbiturates [1:3]
##' length (barbiturates)
##' 
##' barb <- collapse (barbiturates)
##' barb <- orderwl (barb)
##' 
##' plot (barb [1:3], lines.args = list (type = "h"),
##'       col = matlab.dark.palette (3), stacked = TRUE,
##'       stacked.args = list (add.factor = .2))
##' 
##' if (require (latticeExtra)){
##' levelplot (spc ~ .wavelength * z, log (barb), panel = panel.levelplot.points,
##'    cex = 0.3, col = "#00000000", col.regions = matlab.palette (20))
##' }
##' 
##' plotc (apply (barb [,, 42.9~43.2], 1, sum, na.rm = TRUE), spc ~ z,
##'        panel = panel.lines, ylab = expression (I[m/z == 43] / "a.u."))
##' 
NULL
