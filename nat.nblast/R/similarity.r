#' Display two neurons with segments coloured by similarity
#'
#' @param n1 a neuron to compare and colour.
#' @param n2 the neuron to compare against.
#' @param smat a score matrix.
#' @param cols the function to use to colour the segments (e.g. \code{\link{heat.colors}}).
#' @param col the colour with which to draw the comparison neuron.
#' @param AbsoluteScale logical indicating whether the colours should be calculated based on the minimum and maximum similarities for the neuron (\code{AbsoluteScale = FALSE}) or on the minimum and maximum possible for all neurons.
#' @param PlotVectors logical indicating whether the vectors of the \code{dotprops} representation should be plotted. If \code{FALSE}, only the points are plotted.
#' @param ... extra arguments to pass to \code{\link[rgl]{plot3d}}.
#' @return \code{show_similarity} is called for the side effect of drawing the plot; a vector of object IDs is returned.
#' @export
#' @importFrom rgl plot3d
#' @importFrom grDevices colorRampPalette
show_similarity <- function(n1, n2, smat=get(getOption('nat.nblast.defaultsmat')), cols=colorRampPalette(c('#0000FF', '#FF0000')), col='black', AbsoluteScale=FALSE, PlotVectors=TRUE, ...) {
  res <- WeightedNNBasedLinesetMatching.dotprops(n1, n2, NNDistFun=lodsby2dhist, smat=smat, Return='elements')

  if(AbsoluteScale) {
    smat.unique.ordered <- unique(smat)[order(unique(smat))]
    coltable <- rev(cols(length(smat.unique.ordered)))
    segcols <- coltable[sapply(res, function(x) which(x == smat.unique.ordered))]
  } else {
    res.unique.ordered <- unique(res)[order(unique(res))]
    coltable <- rev(cols(length(res.unique.ordered)))
    segcols <- coltable[sapply(res, function(x) which(x == res.unique.ordered))]
  }

  if(PlotVectors) {
    # We need to duplicate each colour as we are drawing line segments, not points
    segcols <- c(sapply(segcols, function(x) c(x,x)))
    plot3d(n2, col=col, PlotVectors=TRUE, PlotPoints=FALSE, ...)
    plot3d(n1, col=segcols, PlotVectors=TRUE, PlotPoints=FALSE, ...)
  } else {
    plot3d(n2, col=col, PlotVectors=FALSE, PlotPoints=TRUE, ...)
    plot3d(n1, col=segcols, PlotVectors=FALSE, PlotPoints=TRUE, ...)
  }
}
