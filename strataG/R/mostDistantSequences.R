#' @title Most Distant Sequences
#' @description Finds the set of sequences that have the greatest mean 
#'   pairwise distance and smallest variance of pairwise distances.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param num.seqs number of sequences to return.
#' @param model a character string specifying the evolutionary model to be used. 
#'   See \link{dist.dna} for more information.
#' @param pairwise.deletion a logical indicating whether to delete sites 
#'   with missing data. See \link{dist.dna} for more information.
#' 
#' @return a vector of the sequence names that are have the greatest mean 
#'   pairwise distance and smallest variance of pairwise distances. 
#'   The names are returned in order from most to least distant.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#' 
#' @examples
#' library(ape)
#' data(dolph.haps)
#' 
#' mostDistantSequences(as.DNAbin(dolph.haps))
#' 
#' @importFrom stats cmdscale var
#' @export
#' 
mostDistantSequences <- function(x, num.seqs = NULL, model = "raw", 
                                 pairwise.deletion = TRUE) {  
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  x <- as.matrix(x)
  if(is.null(num.seqs)) num.seqs <- nrow(x)
  if(num.seqs > nrow(x)) stop(paste("'num.seqs' must be <=", nrow(x)))
  
  # calculate distance between sequences
  seq.dist <- dist.dna(
    x, model = model, pairwise.deletion = pairwise.deletion, as.matrix = TRUE
  )
  
  # convert distances to coordinates
  opt <- options(warn = -1)
  seq.cmd <- cmdscale(seq.dist, k = nrow(x) - 1)
  options(opt)
  
  # normalize coordinates to have mean of 0
  seq.range.mean <- apply(seq.cmd, 2, function(p) mean(range(p, na.rm = TRUE)))
  seq.cmd <- t(t(seq.cmd) - seq.range.mean)
  
  # calculate euclidean distance to mean
  euc.func <- function(p) sqrt(sum(p ^ 2, na.rm = TRUE))
  euc.dist <- sort(apply(seq.cmd, 1, euc.func), decreasing = TRUE)
  
  # get the most distant 'num.seqs' sequences
  ids <- rep(as.character(NA), length = num.seqs)
  ids[1] <- names(euc.dist)[1]
  if(num.seqs == 1) return(ids)
  for(i in 2:num.seqs) {      
    # add sequence with greatest mean distance and smallest variance
    current.ids <- ids[!is.na(ids)]
    other.seqs <- setdiff(rownames(seq.dist), current.ids)
    mean.var <- data.frame(t(sapply(other.seqs, function(hap) {
      m <- mean(seq.dist[hap, current.ids])
      v <- var(seq.dist[hap, current.ids])
      c(mean = m, var = v)
    })))
    mean.var <- cbind(mean.var, euc.dist = euc.dist[rownames(mean.var)])
    ord <- order(mean.var$mean, -mean.var$var, mean.var$euc.dist, 
                 decreasing = TRUE)
    mean.var <- mean.var[ord, ]
    ids[i] <- rownames(mean.var)[1]
  }
  ids
}