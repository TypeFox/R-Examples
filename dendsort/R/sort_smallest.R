#' Sorting and reordering dendrogram nodes by the smallest value
#' 
#' \code{sort_smallest} sorts a dendrogram object based on 
#' the smallest distance in its subtrees, recursively. 
#' The cluster with the smallest distance is placed on the left 
#' side of branch.When a leaf merge with a cluster, the leaf is 
#' placed on the right side.
#'
#' @param d A dendrogram object.\code{d}
#' @return output A sorted dendrogram object. 
#' @keywords internal 
#' @export sort_smallest
#' @aliases sort_smallest
#'  
#' @examples
#' #generate sample data
#' set.seed(1234); par(mar=c(0,0,0,0))
#' x <- rnorm(10, mean=rep(1:5, each=2), sd=0.4)
#' y <- rnorm(10, mean=rep(c(1,2), each=5), sd=0.4)
#' dataFrame <- data.frame(x=x, y=y, row.names=c(1:10))
#' #calculate Euclidian distance
#' distxy <- dist(dataFrame)
#' #hierachical clustering "complete" linkage by default
#' hc <- hclust(distxy)
#' 
#' #sort dendrogram
#' dd <- dendsort(as.dendrogram(hc))
#' hc_sorted  <- as.hclust(dd)
#' 
#' #sort in reverse, you can also pass hclust object
#' plot(dendsort(hc, isReverse=TRUE))
#' 
#' #sort by average distance
#' plot(dendsort(hc, type="average"))
#' 
#' #plot the result
#' par(mfrow = c(1, 3), mai=c(0.8,0.8,2,0.8))
#' plot(x, y, col="gray", pch=19, cex=2)
#' text(x, y, labels=as.character(1:10), cex=0.9)
#' plot(hc,main="before sorting", xlab="", sub="")
#' plot(hc_sorted, main="after sorting", xlab="", sub="")
#' 

sort_smallest <- function(d) {
  if(class(d)!="dendrogram"){stop("d variable must be a dendrogram")}
  
  a <- attributes(d)
  left <- d[[1]]
  right <- d[[2]]
  
  if (is.leaf(left) && is.leaf(right)) {
    n <- d
    attr(n,"min") <- a$height
  } else if (!is.leaf(left) && is.leaf(right)) {
    n <- merge(sort_smallest(left),right,height=a$height)
    attr(n,"min") <- min(a$height, attributes(n[[1]])$min)
  } else if (is.leaf(left) && !is.leaf(right)) {
    n <- merge(sort_smallest(right), left, height=a$height)
    attr(n,"min") <- min(a$height, attributes(n[[1]])$min)
  } else {
    lft  <- sort_smallest(left)
    rght <- sort_smallest(right)
    left.min <- attributes(lft)$min
    right.min <- attributes(rght)$min
    if (left.min <= right.min) {
      n <- merge(lft, rght, height=a$height)
    } else {
      n <- merge(rght, lft, height=a$height)
    }
    attr(n,"min") <- min(a$height, left.min, right.min)
  }
  return(n)
}
