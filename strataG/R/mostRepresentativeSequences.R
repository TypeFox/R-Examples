#' @title Representative Sequences
#' @description Finds the set of sequences that represent the requested 
#'   number of clusters.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param num.seqs number of sequence names to return.
#' @param model a character string specifying the evolutionary model to be used. 
#'   See \link{dist.dna} for more information.
#' @param pairwise.deletion a logical indicating whether to delete sites 
#'   with missing data. See \link{dist.dna} for more information.
#' 
#' @return a vector of the sequence names.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' library(ape)
#' data(dolph.seqs)
#' 
#' mostRepresentativeSequences(as.DNAbin(dolph.seqs))
#' 
#' @importFrom stats cmdscale kmeans
#' @export
#' 
mostRepresentativeSequences <- function(x, num.seqs = NULL, model = "raw", 
                                        pairwise.deletion = TRUE) { 
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  x <- as.matrix(x)
  if(is.null(num.seqs)) num.seqs <- nrow(x)
  num.seqs <- min(num.seqs, nrow(x))
  num.seqs <- max(1, num.seqs)
  if(num.seqs == nrow(x)) return(rownames(x))
  
  # calculate distance between sequences
  seq.dist <- dist.dna(
    x, model = model, pairwise.deletion = pairwise.deletion, as.matrix = TRUE
  )
  
  # convert distances to coordinates
  opt <- options(warn = -1)
  seq.cmd <- cmdscale(seq.dist, k = nrow(x) - 1)
  
  seq.cl <- kmeans(seq.cmd, centers = num.seqs)$cluster
  
  result <- unname(tapply(names(seq.cl), seq.cl, function(ind) {
    if(length(ind) == 1) return(ind)
    if(length(ind) == 2) return(sample(ind, 1))
    ind.dist <- as.matrix(seq.dist[ind, ind])
    cl.cmd <- cmdscale(ind.dist, k = nrow(ind.dist) - 1)
    cmd.mean <- colMeans(cl.cmd)
    dist.to.centroid <- sapply(1:nrow(cl.cmd), function(i) {
      sqrt(sum((cl.cmd[i, ] - cmd.mean) ^ 2))
    })
    ind[which.min(dist.to.centroid)]
  }))
  
  options(opt)
  result[order(match(result, rownames(x)))]
}