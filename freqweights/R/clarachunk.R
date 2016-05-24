## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-09-01 12:36 emilio on emilio-despacho>
## ============================================================




##' Clustering Large Chunks
##'
##' Clustering data splitted in several chunks into k clusters.
##'
##' See \code{\link[cluster]{clara}} for further details.
##'
##' See Examples.
##' @param x data matrix or data frame, each row corresponds to an observation, and each column corresponds to a variable. All variables must be numeric. Missing values (NAs) are allowed.
##' @param k integer, the number of clusters. It is required that 0 < k < n where n is the number of observations of each chunk (i.e., n = nrow(x)).
##' @param samples integer, number of samples to be drawn from the dataset.
##' @return A list with the following values (see \code{\link[cluster]{clara}}):
##' \item{n}{number of rows of the data set.}
##' \item{sample}{labels or case numbers of the observations in the best sample, that is, the sample used by the clara algorithm for the final partition.}
##' \item{medoids}{the medoids or representative objects of the clusters. It is a matrix with in each row the coordinates of one medoid.}
##' \item{tablefreq}{a table of frequency. It is an approximation to the number of cases in each group.}
##' @note This function is based on \code{\link[cluster]{clara}} that it is not available on Windows. Therefore, this implementation does not run on Windows.
##' @seealso \code{\link[cluster]{clara}}, \code{\link{make.readchunk}}
##' @references
##' Antonio Piccolboni \code{mclust.mr} \url{https://github.com/RevolutionAnalytics/rmr2/blob/master/pkg/examples/mclust.mr.R}
##' @examples
##' if(require(cluster)){
##'   k <- 3
##' 
##'   chunk1 <- iris[1:30,1:4]
##'   clus1 <- clarasub(chunk1,k)
##' 
##'   chunk2 <- iris[-c(1:30),1:4]
##'   clus2 <- clarasub(chunk2,k)
##' 
##'   subclusters <- list(clus1, clus2)
##'   b <- claramerge(subclusters,k)
##'      print(b$medoids)
##'
##'    print(nrow(b$tablefreq))
##'   print(b$tablefreq)
##' }
##' @name clarachunk
##' @rdname clarachunk
NULL


##' @rdname clarachunk
##' @export
clarasub <- function(x, k, samples = 50){
  ## CPU time O( n * p * samplesize * samples). Storage O( n * p) + O( samplesize)
##   if(require("cluster")){
##   clus <- clara(x, k, samples = samples, rngR=TRUE, keep.data = FALSE)
##   return(list(n=NROW(x), sample=x[clus$sample,,drop=FALSE], medoids= clus$medoids, tablefreq=tablefreq(clus$clustering)))
## }
  if (requireNamespace("cluster", quietly = TRUE)) {
    clus <- cluster::clara(x, k, samples = samples, rngR=TRUE, keep.data = FALSE)
    return(list(n=NROW(x), sample=x[clus$sample,,drop=FALSE], medoids= clus$medoids, tablefreq=tablefreq(clus$clustering)))
  } else {
      ## do something else not involving rgl.
   }
  
}


##' @param subclusters list of objects returned by clarasub 
##' @rdname clarachunk
##' @export
claramerge <- function(subclusters, k, samples= 50){
  if(!is.null(names(subclusters)) && all(names(subclusters) == c("n","sample","medoids","tablefreq"))) {
    ## just one cluster
    return(subclusters)
  }
  sizes <- laply(subclusters, function(x)x[["n"]])
  totalsize <- sum(sizes)
  rangesize <- range(sizes)
  ratiosize <- max(rangesize)/min(rangesize)
  resample <- function(x) x$sample[ sample(seq_len(NROW(x$sample)), round( NROW(x$sample) * ratiosize), replace=TRUE),,drop=FALSE]
  clus <- clarasub( x=ldply(subclusters, resample), k=k, samples = samples)
  clus$n <- totalsize
  clus$tablefreq <-  tablefreq(ldply(subclusters, function(x) x$tablefreq), freq="freq")
  clus
}

