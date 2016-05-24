#' Calculate the x coordinates given a branch of dendrogram
#' 
#' \code{cal_node_limit} is a code modified from plotNodeLimit()
#' to x coordinates of branches given a branch of dendrogram.
#' 
#' @param x1 A x coordinate\code{x1}
#' @param x2 Another x coordinate\code{x2}
#' @param subtree A dendrogram object.\code{subtree}
#' @param center A logical whether the dendrogram is centered.\code{center}
#' 
#' @return output A list of parameters. 
#'
#' @keywords internal 
#'
#' @export cal_node_limit
#' @aliases cal_node_limit
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
#' total <- cal_total_length(as.dendrogram(hc))
#' 

cal_node_limit <- function(x1, x2, subtree, center){
  ## get the left borders limit[k] of all children k=1..K, and
  ## the handle point `x' for the edge connecting to the parent.
  inner <- !is.leaf(subtree) && x1 != x2
  if(inner) {
    K <- length(subtree)
    mTop <- .memberDend(subtree)
    limit <- integer(K)
    xx1 <- x1
    for(k in 1L:K) {
      m <- .memberDend(subtree[[k]])
      ##if(is.null(m)) m <- 1
      xx1 <- xx1 + (if(center) (x2-x1) * m/mTop else m)
      limit[k] <- xx1
    }
    limit <- c(x1, limit)
  } else { ## leaf
    limit <- c(x1, x2)
  }
  mid <- attr(subtree, "midpoint")
  center <- center || (inner && !is.numeric(mid))
  x <- if(center) mean(c(x1,x2)) else x1 + (if(inner) mid else 0)
  list(x = x, limit = limit)
}