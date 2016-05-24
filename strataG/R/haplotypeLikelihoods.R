#' @title Haplotype Likelihoods
#' @description Calculate likelihood of each haplotype based on gamma 
#'   distribution of pairwise distances.
#' 
#' @param x a \code{\link[ape]{DNAbin}} object.
#' @param model a character string specifying the evolutionary model to be 
#'   used. Passed to \code{\link[ape]{dist.dna}}. 
#' @param pairwise.deletion a logical indicating whether to delete the 
#'   sites with missing data in a pairwise way. Passed to 
#'   \code{\link[ape]{dist.dna}}.
#' @param plot.n number of haplotypes with lowest delta(log-likelihoods) to 
#'   plot. Defaults to all haplotypes. Set to 0 to supress plotting.
#' @param ... arguments passed from other functions (ignored).
#' 
#' @details Fits a Gamma distribution to the pairwise distances of haplotypes 
#'   and calculates the log-likelihood for each haplotype (sum of all pairwise 
#'   log-likelihoods for that haplotype). Haplotypes that are extremely 
#'   different from all others will have low log-likelihoods. Values returned 
#'   as delta(log-likelhoods) = difference of log-likelihoods from maximum 
#'   observed values.
#' 
#' @return vector of delta(log-Likelihoods) for each haplotype, sorted from 
#'   smallest to largest, and a plot of their distributions.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov}
#'
#' @examples
#' library(ape)
#' data(dolph.haps)
#' 
#' haplotypeLikelihoods(as.DNAbin(dolph.haps))
#' 
#' @importFrom stats sd dgamma
#' @importFrom graphics dotchart
#' @export
#' 
haplotypeLikelihoods <- function(x, model = "raw", pairwise.deletion = FALSE, 
                                 plot.n = NULL, ...) {
  
  if(!inherits(x, "DNAbin")) stop("'x' must be a DNAbin object")
  
  # calculate distance between sequences
  seq.dist <- dist.dna(
    x, model = model, pairwise.deletion = pairwise.deletion, as.matrix = TRUE
  )
  if(any(is.nan(seq.dist)) | any(is.na(seq.dist))) {
    warning("NA/NaN in pairwise distance matrix. NULL returned.")
    return(NULL)
  }
  
  dist.vec <- seq.dist[lower.tri(seq.dist)]
  dist.mean <- mean(dist.vec)
  dist.sd <- sd(dist.vec)
  scale <- (dist.sd ^ 2) / dist.mean
  shape <- (dist.mean / dist.sd) ^ 2
  
  log.lik <- sapply(rownames(seq.dist), function(this.seq) {
    this.dist <- seq.dist[this.seq, ]
    this.dist <- this.dist[names(this.dist) != this.seq]
    sum(log(dgamma(this.dist, shape = shape, scale = scale)), na.rm = TRUE)
  })
  delta.log.lik <- sort(log.lik - max(log.lik, na.rm = T), decreasing = FALSE)
  
  if(is.null(plot.n)) plot.n <- length(delta.log.lik)
  n <- min(c(plot.n, length(delta.log.lik)))
  if(plot.n > 0) {
    dotchart(rev(delta.log.lik[1:n]), pch = 19, bg = "black", 
             xlab = expression(paste(Delta, "lnL")), 
             main = "Haplotype delta(log-Likelihoods)"
    )
  }
  
  delta.log.lik
}