#' Ensure an edgelist has all dyads and a column of weights.
#'
#' Given a matrix or data.frame edgelist, fill in all possible edges not already listed with a weight of 0 or the value of \code{fillBlanksWith}.
#' 
#' @param elist data.frame or matrix, see 'Details' for formatting assumptions.
#' @param fillBlanksWith numeric, default weight for edges not already listed in elist.
#' @param nodelist character, optional list of node names.
#' @return data.frame, full list of all possible edges with weights for each in third column.
#' @export
#' @references
#' \url{https://github.com/shaptonstahl/}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @details The \code{elist} can be either a data.frame or a matrix with either 2 or 3 columns. Each row is an edge. The first column lists the node the edge is 'from' and the second column lists the node the edge is 'to'. If there is a third column, it lists the weight of the edge.
#' @examples
#' g <- erdos.renyi.game(10, 2/10)
#' EdgelistFill(get.edgelist(g))
#' EdgelistFill(get.edgelist(g), nodelist=1:10)
#' 
#' E(g)$weight <- runif(ecount(g))
#' el <- cbind(get.edgelist(g), E(g)$weight)
#' EdgelistFill(el)
#' EdgelistFill(el, nodelist=1:10)
EdgelistFill <- function(elist, 
                         fillBlanksWith=0, 
                         nodelist) {
  # Guardians
  if( (is(elist, "matrix") || is(elist, "data.frame")) && 2 == ncol(elist) ) {
    # elist has no weights; perhaps it was the result of get.edgelist(g)
    elist <- cbind(elist, 1)  # put 1s in for each edge
  } else if( (is(elist, "matrix") || is(elist, "data.frame")) && 3 == ncol(elist) ) {
    # elist has a third column, presumably with weights
    stopifnot( is.numeric(elist[,3]) )
  } else {
    stop("elist must have two or three columns")
  }
  
  if( missing(nodelist) ) {
    nodelist <- sort(union(unique(elist[,1]), unique(elist[,2])))
  } else {
    stopifnot(all(elist[,1] %in% nodelist),
              all(elist[,2] %in% nodelist))
    nodelist <- sort(nodelist)
  }
  
  out <- expand.grid(nodelist, nodelist)
  out <- data.frame(out[,2], out[,1], fillBlanksWith)  # ensures sorted by first column, then second column
  if( is.null(names(elist)) ) {
    names(out) <- c("fromnode", "tonode", "weight")
  } else {
    names(out) <- names(elist)
  }
  
  for(i in 1:nrow(elist)) {
    if( is.factor(elist[,1]) ) {
      match1 <- levels(elist[i,1])[as.numeric(elist[i,1])]  # "unfactor" in case factor levels don't agree
    } else {
      match1 <- elist[i,1]
    }
    if( is.factor(elist[,2]) ) {
      match2 <- levels(elist[i,2])[as.numeric(elist[i,2])]  # "unfactor" in case factor levels don't agree
    } else {
      match2 <- elist[i,2]
    }
    out[out[,1]==match1 & out[,2]==match2, 3] <- elist[i,3]
  }
  return( out )
}
