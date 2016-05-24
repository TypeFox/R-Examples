#' Plot method for class 'lsbclust'
#' 
#' This plot method simply plots each of the components in the list of class \code{lsbclust}. 
#' 
#' @param x An object of class \code{orc.kmeans}
#' @param type A character vector indicating which component(s) of \code{x} to plot: a combination of
#' \code{"overall"}, \code{"rows"}, \code{"columns"} and \code{"interactions"}.
#' @param biplot.axes A logical indicating whether to plot calibrated biplot axes for the line
#' segments indicated in \code{segments} or not.
#' @param \dots additional arguments passed to the plot methods of the respective components, typically
#' to \code{\link{theme}}. Use e.g. \code{plot(x$interactions)} for more control over the respective
#' plots.
#' @keywords hplot
#' @author Pieter C. Schoonees
#' @method plot lsbclust
#' @seealso \code{\link{plot.int.lsbclust}}, \code{\link{plot.ovl.kmeans}}, 
#' \code{\link{plot.row.kmeans}}, \code{\link{plot.col.kmeans}}
#' @export
#' @examples
#' data("dcars")
#' m <- lsbclust(data = dcars, margin = 3, delta = c(1, 1, 1, 1), nclust = 5, nstart = 1)
#' plot(m)
plot.lsbclust <- function(x, type = c("overall", "rows", "columns", "interactions"), 
                          biplot.axes = TRUE, ...) {
  
  ## Check which type of plots to produce
  type <- match.arg(tolower(type), choices = c("overall", "rows", "columns", "interactions"), several.ok = TRUE)
  type <- type[type %in% names(x)[!sapply(x, is.null)]]
  
  ## Setup list of grobs
  plots <- vector(mode = "list", length = 4)
  names(plots) <- c("overall", "rows", "columns", "interactions")
  
  ## Fill in the list
  if ("overall" %in% type) {
    plots[["overall"]] <- plot(x[["overall"]], ...)
  }
  if ("rows" %in% type) {
    plots[["rows"]] <- plot(x[["rows"]], ...)
  }
  if ("columns" %in% type) {
    plots[["columns"]] <- plot(x[["columns"]], ...)
  }
  if ("interactions" %in% type) {
    plots[["interactions"]] <- plot(x[["interactions"]], biplot.axes = biplot.axes, ...)
  }
  return(plots[type])
}