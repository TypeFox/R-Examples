#######################################################################
##
## Function: plot.anchors.rank()
## Author  : Jonathan Wand <wand@stanford.edu>
##
#######################################################################

plot.anchors.rank <- function(x, ... ,
                              xy   = c("minimum","interval")) {

  if (is.null(x$combn))
    stop("plot.anchors.anchors is currently available only if anchors(...,combn=TRUE) was used;",
         "perhaps you are looking for 'barplot.anchors.rank' to plot the disribution of ranks")

  plot.anchors.combn( x$combn, ..., xy=xy)
  
}
