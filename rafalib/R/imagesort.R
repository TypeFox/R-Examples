#' image with sorted rows
#'
#' the rows are sorted such that the first column has 2 blocks,
#' the second column has 4 blocks, etc. see example("imagesort")
#' 
#' @param x a matrix of 0s and 1s
#' @param col the colors of 0 and 1
#' @param ... arguments to heatmap
#' 
#' @author Michael I. Love
#' 
#' @examples
#'
#' x <- replicate(4,sample(0:1,40,TRUE))
#' imagesort(x)
#' 
imagesort <- function(x,col=c("white","black"),...) {
  hc <- hclust(dist(t(x)))
  y <- sweep(x,2,2^(ncol(x)-order(hc$order)),"*")
  z <- x[order(rowSums(y)),]
  heatmap(z, Rowv=NA, 
          Colv=as.dendrogram(hc),
          labRow=FALSE,
          scale="none",
          col=col,...)
}
