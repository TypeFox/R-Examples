##' Color coding legend for factors
##' Modifies a list of lattice arguments (as for \code{\link[lattice]{levelplot}}, etc.) according to
##' the factor levels. The colorkey will shows all levels (including unused), and the drawing colors
##' will be set accordingly.
##' 
##' \code{trellis.factor.key} is used during \code{levelplot}-based plotting of factors (for
##' hyperSpec objects) unless \code{transform.factor = FALSE} is specified.
##' 
##' @param f the factor that will be color-coded
##' @param levelplot.args a list with levelplot arguments
##' @return the modified list with levelplot arguments.
##' @author C. Beleites
##' @seealso \code{\link[lattice]{levelplot}}
##' @keywords aplot
##' @export
##' @importFrom lattice level.colors
##' @examples
##' 
##' chondro$z <- factor (rep (c("a", "a", "d", "c"),
##'                           length.out = nrow (chondro)),
##'                      levels = letters [1 : 4])
##' 
##' str (trellis.factor.key (chondro$z))
##' 
##' plotmap (chondro, z ~ x * y)
##' 
##' ## switch off using trellis.factor.key:
##' ## note that the factor levels are collapsed to c(1, 2, 3) rather than
##' ## c (1, 3, 4)
##' plotmap (chondro, z ~ x * y, transform.factor = FALSE)
##' 
##' plotmap (chondro, z ~ x * y,
##'          col.regions = c ("gray", "red", "blue", "dark green"))
##' 
trellis.factor.key <- function (f, levelplot.args = list ()) {
  at <-  seq (0, nlevels (f)) + .5
  
  if (is.null (levelplot.args$col.regions))
    cols <- level.colors (seq_along (levels (f)), at)
  else
    cols <- level.colors (seq_along (levels (f)), at, levelplot.args$col.regions)
    
  modifyList (list (at = at,
                    col.regions = cols, 
                    colorkey = list (lab = list (at = seq_along (levels (f)),
                                                 lab = levels (f)))),
              levelplot.args)
  
}
