##' This function implements a version of the hierarchical, agglomerative clustering  \code{\link[fastcluster]{hclust.vector}} focused on table of frequencies.
##' 
##' Any variables in the formula are removed from the data set.
##' 
##' This function is a wrapper of \code{\link[fastcluster]{hclust.vector}} to be used with tables of frequencies. It use the frequency weights as parameter \code{members}.
##' @title Fast hierarchical, agglomerative clustering of frequency data
##' @param data any object that can be coerced into a double matrix
##' @param method the agglomeration method to be used. This must be (an unambiguous abbreviation of) one of "\code{single}", "\code{ward}", "\code{centroid}" or "\code{median}".
##' @param freq a one-sided, single term formula specifying frequency weights
##' @param metric the distance measure to be used. This must be one of \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} or \code{"minkowski"}
##' @param p parameter for the Minkowski metric.
##' @seealso \code{\link[fastcluster]{hclust.vector}}, \code{link{tablefreq}}
##' @importFrom fastcluster hclust hclust.vector
##' @import dplyr
##' @export
##' @rdname hclustvfreq
##' @examples
##' library(dplyr)
##' library(fastcluster)
##' 
##' data <- iris[,1:3,drop=FALSE]
##' hc <- hclustvfreq(data, method="centroid",metric="euclidean")
##' cutree(hc,3) ## Different length than data
##'
##' tfq <- tablefreq(iris[,1:3])
##' hc <- .hclustvfreq(tfq, method="centroid",metric="euclidean")
##' tfq$group <- cutree(hc,3)
hclustvfreq <- function(data, freq=NULL, method="single", metric= "euclidean", p=NULL){
   tablefreq(data,freq=freq) %>%
     .hclustvfreq(method=method, metric=metric, p=p)
}


##' @param tfq a frequency table
##' @rdname hclustvfreq
##' @export
.hclustvfreq<- function(tfq,method="single", metric= "euclidean", p=NULL){
  if(!all(complete.cases(tfq)) ) {
    stop("hclustvfreq: Missing cases, or non-positive frequency weights!")
  }
  hclust.vector(tfq[,-ncol(tfq),drop=FALSE], method=method, members=unlist(tfq[,ncol(tfq)]), metric=metric, p=p)
}
