#' plot clusters
#' @param x spatcluster-cluster object
#' @param data point pattern object used for computing the graph
#' @param add add or plot new
#' @param atleast plot only cluster with 'atleast' points in them
#' @param col colors for clusters, chosen randomly if missing.
#' @param ... passed to points
#'
#'
#' @export

plot.sgc<-function(x, data, atleast=2, add=FALSE, col, ...)
{
  data <- sg_parse_coordinates(data)
  w <- (1:x$nclusters)[sapply(x$clusters,length)>=atleast]
  n <- length(w)

  if(missing(col)) col <- rgb(red=runif(n,0.1,1),green=runif(n,0.1,1), blue=runif(n,0.1,1) )
  if(ncol(data)==2) {
      if(!add) plot(NA, NA, xlim=range(data[,1]), ylim=range(data[,2]), asp=1, xlab="x", ylab="y")
      for(i in w)
      {
        points(data[x$clusters[[i]],1], data[x$clusters[[i]], 2], col=col[i], ...)
      }
  }
  else stop("Plotting of clusters only in 2D")

}
