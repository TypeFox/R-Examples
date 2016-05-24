#################################################################################
###
###  plotvoronoi - plot spectral maps with irregular point pattern
###  
###  plots intensity or extra data column over 2 extra data columns

##' @param use.tripack Whether package tripack should be used for calculating
##'   the voronoi polygons. If \code{FALSE}, package deldir is used instead.
##'   See details.
##' @param mix For Voronoi plots using package tripack, I experienced errors if
##'   the data was spatially ordered. Randomly rearrangig the rows of the
##'   hyperSpec object circumvents this problem.
##' @rdname levelplot
##' @include levelplot.R
##' @export
##' @seealso \code{\link[latticeExtra]{panel.voronoi}}
##' @importFrom latticeExtra panel.voronoi
##' @importFrom lattice prepanel.default.levelplot
##' 
plotvoronoi <- function (object, model = spc ~ x * y,
                         use.tripack = FALSE, mix = FALSE, ...){
  if (!requireNamespace ("latticeExtra"))
   stop ("package latticeExtra is needed for Voronoi plots.")

  if (use.tripack){
    if (!requireNamespace ("tripack"))
      stop ("package tripack requested but not available.")
  } else {
    if (!requireNamespace ("deldir"))
      stop ("package deldir requested but not available.")
  }
 
  if (use.tripack && mix)
      object@data <- object@data [sample (nrow (object)),]

  dots <- modifyList (list (object = object,
                            model = model,
                            panel = panel.voronoi,
                            prepanel = prepanel.default.levelplot,
                            pch = 19, cex = .25,
                            col.symbol = "#00000020",
                            border = "#00000020",
                            use.tripack = use.tripack),
                      list (...))
  do.call (plotmap, dots)
}
