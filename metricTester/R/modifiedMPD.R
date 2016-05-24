#' Calculate different versions of abundance-weighted MPD
#'
#' Given a picante-style community data matrix (sites are rows, species are columns), 
#' a phylogenetic distance matrix, and a desired method of abundance-weighting, will
#' calculate MPD.
#'
#' @param samp A picante-style community data matrix with sites as rows, and
#' species as columns
#' @param dis Phylogenetic distance matrix
#' @param abundance.weighted One of either "FALSE", "interspecific",
#' "intraspecific", or "complete"
#' 
#' @details See accompanying publication for details. Non-abundance-weighted and
#' interspecific and intraspecific methods are equivalent to those previously described by
#' Clarke & Warwick.
#'
#' @return A vector of MPD values, calculated according to the abundance-weighted method
#' specified.
#'
#' @export
#'
#' @importFrom ape cophenetic.phylo
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' #simulate tree with birth-death process
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' dists <- ape::cophenetic.phylo(tree)
#'
#' results <- modifiedMPD(cdm, dists, abundance.weighted = "interspecific")

modifiedMPD <- function(samp, dis, abundance.weighted = FALSE) 
{
    N <- dim(samp)[1]
    mpd <- numeric(N)
    for (i in 1:N) {
        sppInSample <- names(samp[i, samp[i, ] > 0])
        if (length(sppInSample) > 1) {
            sample.dis <- dis[sppInSample, sppInSample]
            if (abundance.weighted == "interspecific") {
                sample.weights <- t(as.matrix(samp[i, sppInSample, 
                  drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                  drop = FALSE])
	            diag(sample.weights) <- 0
                mpd[i] <- weighted.mean(sample.dis, sample.weights)
            }
            else if (abundance.weighted == "intraspecific") {
                sample.weights <- t(as.matrix(samp[i, sppInSample, 
                  drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                  drop = FALSE])
	            diag(sample.weights) <- diag(sample.weights) - sqrt(diag(sample.weights))
                mpd[i] <- weighted.mean(sample.dis, sample.weights)
            }
            else if (abundance.weighted == "complete") {
                sample.weights <- t(as.matrix(samp[i, sppInSample, 
                  drop = FALSE])) %*% as.matrix(samp[i, sppInSample, 
                  drop = FALSE])
                mpd[i] <- weighted.mean(sample.dis, sample.weights)
            }
            else {
                mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
            }
        }
        else {
            mpd[i] <- NA
        }
    }
    mpd
}
