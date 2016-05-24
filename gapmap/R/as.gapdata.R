#'Make a gapdata class object
#'
#'This function generates a gapdata class object.  This object is used for drawing dendrograms and heatmaps.
#'
#' @param d dendrogram class object
#' @param segments a data.frame containing segments information
#' @param labels a data.frame containing labels information
#' @param labels_df data.frame storking the label positions
#' @param ... ignored
#' @export as.gapdata
#' @aliases as.gapdata
#' @return the gapdata class object
#' @keywords internal
#' 

#make a gapdata class
as.gapdata <- function(d, segments, labels, ...){
  x <- list(dendrogram = d, segments = segments, labels=labels)
  class(x) <- "gapdata"
  x
}