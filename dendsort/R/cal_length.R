#' Recursive function to calculate the length of branches
#' 
#' \code{cal_length} is a code modified from plotNode()
#' to calculate the length of lines to draw the branch of a dendrogram. This
#' function was developed to evaluate the use of ink for visualization.
#' 
#' @param x1 A x coordinate\code{x1}
#' @param x2 Another x coordinate\code{x2}
#' @param subtree A dendrogram object.\code{subtree}
#' @param center A logical whether the dendrogram is centered.\code{center}
#' @param nodePar A node parameter.\code{nodePar}
#' @param edgePar An edge parameter.\code{edgePar}
#' @param horiz A logical about layout.\code{horiz}
#' @param sum A sum of length.\code{sum}
#'
#' @return output The length. 
#'
#' @keywords internal 
#'
#' @export cal_length
#' @aliases cal_length
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
cal_length <-function(x1, x2, subtree, center, nodePar, edgePar, horiz = FALSE, sum){
  inner <- !is.leaf(subtree) && x1 != x2
  yTop <- attr(subtree, "height")
  bx <- cal_node_limit(x1, x2, subtree, center)
  xTop <- bx$x
  
  ## handle node specific parameters in "nodePar":
  hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
  if(!hasP) nPar <- nodePar
  
  if (is.leaf(subtree)) {
  }else if (inner) {
    for (k in seq_along(subtree)) {
      child <- subtree[[k]]
      ## draw lines to the children and draw them recursively
      yBot <- attr(child, "height")
      if (getOption("verbose")) cat("ch.", k, "@ h=", yBot, "; ")
      if (is.null(yBot))
        yBot <- 0
      xBot <-
        if (center) mean(bx$limit[k:(k + 1)])
      else bx$limit[k] + .midDend(child)
      
      hasE <- !is.null(ePar <- attr(child, "edgePar"))
      if (!hasE)
        ePar <- edgePar
      i <- if (!is.leaf(child) || hasE) 1 else 2
      
      #calculate lentth
      len = abs(xTop-xBot) + abs(yTop-yBot)
      #         cat(c("\n", len))
      vln <- NULL
      sum = sum + len; 
      sum = cal_length(bx$limit[k], bx$limit[k + 1], subtree = child, center, nodePar, edgePar, horiz, sum)
    }
  }
  #     invisible()
  return(sum)
}