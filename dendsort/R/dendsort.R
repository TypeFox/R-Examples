#' Sorting and reordering dendrogram nodes
#' 
#' \code{dendsort} sorts a dendrogram object which is 
#' typically  a result of hierarchical clustering (hclust). The 
#' subtrees in the resulting dendrogram are sorted based on the 
#' average distance of subtrees at every merging point. The 
#' tighter cluster, in other words the cluster with smaller 
#' average distance, is placed on the left side of branch.  
#' When a leaf merge with a cluster, the leaf is placed on the 
#' right side.
#'
#' @param d a dendrogram or hclust object.\code{d}
#' @param isReverse logical indicating if the order should be reversed.Defaults to FALSE\code{isReverse}
#' @param type character indicating the type of sorting. Default to "min" \code{type}
#' @return output A sorted dendrogram or hclust. 
#' @keywords dendrogram 
#' @export dendsort
#' @aliases dendsort
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

dendsort <- function(d, isReverse=FALSE, type="min") {
  if(!inherits(d, "dendrogram") && !inherits(d, "hclust")){
    stop("d variable must be a dendrogram or hclust object")
  }
  
  #assign dendrogram
  dend = d
  if(inherits(d, "hclust")){
   dend = as.dendrogram(d)
  }
  #type string to lower case
  type = tolower(type)
  if(type=="average"){
    #sort by average distance
    if(isReverse){
      #sort in reverse
      n = sort_average_r(dend)
    }else{
      #sort in left to right order 
      n = sort_average(dend)
    }
  }else if(type == "min"){
    #sort by smallest distance
    if(isReverse){
      n = sort_smallest_r(dend)
    }else{
      n = sort_smallest(dend)
    }
  }else{
    stop("unrecognized type variable "+type)
  }
  #if input was a hclust object, convert it back to hclust
  if(inherits(d, "hclust")){
    n = as.hclust(n)
  }
  
  return(n)
}

