#' Calculate the total length of lines to draw the dendrogram
#' 
#' \code{cal_total_length} is a code modified from plot.dendrogram()
#' to calculate the total length of lines to draw a dendrogram. This
#' function was developed to evaluate the use of ink for visualization.
#'
#' @param x A dendrogram object.\code{x}
#'
#' @return output The total length. 
#'
#' @keywords internal 
#'
#' @export cal_total_length
#' @aliases cal_total_length
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
#' total_length <- cal_total_length(as.dendrogram(hc))
#' 

#modifided from plot.dendrogram()
cal_total_length <-function (x){
  center = FALSE
  edge.root = is.leaf(x) || !is.null(attr(x, "edgetext"))
  nodePar = NULL
  edgePar = list()
  dLeaf = NULL
  xlab = ""
  ylab = ""
  xaxt="n"
  yaxt="s"
  horiz = FALSE
  
  hgt <- attr(x, "height")
  if (edge.root && is.logical(edge.root))
    edge.root <- 0.0625 * if(is.leaf(x)) 1 else hgt
  mem.x <- .memberDend(x)
  yTop <- hgt + edge.root
  if(center) { x1 <- 0.5 ; x2 <- mem.x + 0.5 }
  else       { x1 <- 1   ; x2 <- mem.x }
  xl. <- c(x1 - 1/2, x2 + 1/2)
  yl. <- c(0, yTop)
  if (horiz) {## swap and reverse direction on `x':
    tmp <- xl.; xl. <- rev(yl.); yl. <- tmp
    tmp <- xaxt; xaxt <- yaxt; yaxt <- tmp
  }
  
  l = cal_length(x1, x2, x, center = center, nodePar = nodePar, edgePar = edgePar, horiz = horiz, sum=0)
  return(l)
}




.memberDend <- function(x) {
  r <- attr(x,"x.member")
  if(is.null(r)) {
    r <- attr(x,"members")
    if(is.null(r)) r <- 1L
  }
  r
}
.midDend <- function(x){
  if(is.null(mp <- attr(x, "midpoint"))) 0 else mp
}


