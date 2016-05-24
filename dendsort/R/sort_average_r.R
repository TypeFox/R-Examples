#' Sorting and reordering dendrogram nodes by average distances in reverse
#' 
#' \code{sort_average_r} sorts a dendrogram object in reverse based on 
#' the average distance of its subtrees, recursively. 
#' The tighter cluster, in other words the cluster with smaller 
#' average distance, is placed on the right side of branch.  
#' When a leaf merge with a cluster, the leaf is placed on the 
#' left side.
#'
#' @param d A dendrogram object.\code{d}
#'
#' @return output A sorted dendrogram object. 
#'
#' @keywords internal 
#'
#' @export sort_average_r
#' @aliases sort_average_r
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
sort_average_r <- function(d) {
  if(class(d)!="dendrogram"){stop("d variable must be a dendrogram")}
  
  a <- attributes(d)
  left <- d[[1]]
  right <- d[[2]]
  
  if (is.leaf(left) && is.leaf(right)) {
    n <- d
    attr(n,"sum") <- a$height
  } else if (!is.leaf(left) && is.leaf(right)) {
    n <- merge(right, sort_average_r(left),height=a$height)
    attr(n,"sum") <- a$height + attributes(n[[2]])$sum
  } else if (is.leaf(left) && !is.leaf(right)) {
    n <- merge(left, sort_average_r(right), height=a$height)
    attr(n,"sum") <- a$height + attributes(n[[2]])$sum
  } else {
    lft  <- sort_average_r(left)
    rght <- sort_average_r(right)
    left.count <- attributes(lft)$members -1
    right.count <- attributes(rght)$members -1
    left.sum <- attributes(lft)$sum
    right.sum <- attributes(rght)$sum
    left.avg <- left.sum / left.count
    right.avg <- right.sum / right.count
    
    if (left.avg >= right.avg) {
      n <- merge(lft, rght, height=a$height)
    } else {
      n <- merge(rght, lft, height=a$height)
    }
    attr(n,"sum") <- a$height + left.sum + right.sum
  }
  return(n)
}
