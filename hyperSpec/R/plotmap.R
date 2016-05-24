#################################################################################
###
###  plotmap - plot spectral maps
###  
###  plots intensity or extra data column over 2 extra data columns

## TODO: check wheter func should be applied or not



##' Plot a Map and Identify/Select Spectra in the Map
##' \code{\link[lattice]{levelplot}} functions for hyperSpec objects.  An image or map of a summary
##' value of each spectrum is plotted. Spectra may be identified by mouse click.
##' 
##' The \code{model} can contain the special column name \code{.wavelength} to specify the wavelength
##' axis.
##' 
##' \code{plotmap}, \code{map.identify}, and the \code{levelplot} methods internally use the same
##' gateway function to \code{\link[lattice]{levelplot}}. Thus \code{transform.factor} can be used
##' with all of them and the panel function defaults to \code{\link[lattice]{panel.levelplot.raster}}
##' for all three. Two special column names, \code{.rownames} and \code{.wavelength} may be used. 
##' 
##' \code{levelplot} plots the spectra matrix.
##' 
##' \code{plotvoronoi} calls \code{plotmap} with different default settings, namely the panel
##' function defaults to \code{\link[latticeExtra]{panel.voronoi}}.
##' \code{\link[latticeExtra]{panel.voronoi}} depends on either of the packages 'tripack' or 'deldir'
##' being installed. For further information, please consult the help page of
##' \code{\link[latticeExtra]{panel.voronoi}}.  On the \code{\link{chondro}} data set, \code{plotmap}
##' is roughly 5 times faster than \code{plotvoronoi} using tripack, and ca. 15 times faster than
##' \code{plotvoronoi} using deldir. Package tripack, however, is free only for non-commercial
##' use. Also, it seems that tripack version hang (R running at full CPU power, but not responding
##' nor finishing the calculation) for certain data sets. In this case, \code{mix = TRUE} may help.
##' 
##' \code{map.identify} calls \code{plotmap} and \code{plotvoronoi}, respectively and waits for
##' (left) mouse clicks on points. Other mouse clicks end the input.
##' 
##' Unlike \code{\link[lattice]{panel.identify}}, the indices returned by \code{map.identify} are in
##' the same order as the points were clicked. Also, multiple clicks on the same point are returned
##' as multiple entries with the same index.
##' 
##' \code{map.identify} uses option \code{debuglevel} similar to \code{\link{spc.identify}}:
##' \code{debuglevel == 1} will plot the tolerance window if no data point was inside (and
##' additionally labels the point) while \code{debuglevel == 2} will always plot the tolerance
##' window.
##' 
##' The \code{map.sel.*} functions offer further interactive selection, see
##' \code{\link{map.sel.poly}}.
##'
##' @rdname levelplot
##' @aliases plotmap plotvoronoi levelplot,formula,hyperSpec-method
##'   levelplot,hyperSpec,missing-method map.identify
##' @param object,data the \code{hyperSpec} object
##' @param model,x formula specifying the columns of object that are to be
##'   displayed by \code{\link[lattice]{levelplot}}
##' @param func,func.args Before plotting, \code{plotmap} applies function
##'   \code{func} with the arguments given in the list \code{func.args} to each
##'   of the spectra. Thus a single summary value is displayed for each of the
##'   spectra.
##'
##' This can be suppressed manually by setting \code{func} to NULL. It is automatically suppressed if
##' \code{.wavelength} appears in the formula.
##' @param voronoi Should the plot for identifying spectra by mouse click be
##'   produced by \code{plotmap} (default) or \code{plotvoronoi}?
##' @param ... further arguments are passed down the call chain, and finally
##'   to \code{\link[lattice]{levelplot}}
##' @return \code{map.identify} returns a vector of row indices into
##'   \code{object} of the clicked points.
##' 
##' The other functions return a lattice object.
##' @author C. Beleites
##' @seealso \code{vignette (plotting)}, \code{vignette (introduction)}
##' 
##' \code{\link{plot}}
##' @export
##' @keywords hplot
##' @examples
##' 
##' \dontrun{
##' vignette (plotting)
##' vignette (introduction)
##' }
##' 
##' levelplot (spc ~ y * x, chondro [,,1003]) # properly rotated
##' plotmap (chondro [,,1003])
##' 
##' # plot spectra matrix
##' levelplot (spc ~ .wavelength * t, laser, contour = TRUE, col = "#00000080")
##' # see also plotmat
##' 
##' plotmap (chondro, clusters ~ x * y)
##' 
##' # Voronoi plots
##' smpl <- sample (chondro, 300)
##' plotmap (smpl, clusters ~ x * y)
##' if (require (tripack)) 
##'     plotvoronoi (smpl, clusters ~ x * y)
##' if (require (deldir)) 
##'     plotvoronoi (smpl, clusters ~ x * y,
##'                  use.tripack = FALSE)
##' 
plotmap <- function (object, model = spc ~ x * y,
                     func = mean, func.args = list (), ...){
  chk.hy (object)
  validObject (object)

  if (! is.null (func) & ! any (grepl ("[.]wavelength", model)))
    object <- do.call (apply, c (list (object, 1, func), func.args))
  
  dots <- modifyList (list (aspect = "iso"),
                      list (...))
                       
  dots <- c (list (x = model, data = object), dots)

  do.call(.levelplot, dots)
}

