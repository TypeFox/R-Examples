#' Pathway Level Statistics
#' 
#' Takes a vector of statistics with each element corresponds to a gene or phosphorylation site, 
#' and calculates a combined statistics for those that belong to the same pathway or kinase.
#' 
#' @usage pathwayStats(PGs, T, minSize=5, method="Stouffer")
#' @param PGs an array of names indicating genes or substrates that belong to a given pathway or kinase.
#' @param T a vector of statistics (z-scores converted) with each element correspond to a gene or 
#' phosphorylation site that belong to the same pathway or kinase.
#' @param minSize the size of annotation groups to be considered for calculating enrichment. Groups 
#' that are smaller than the minSize will be removed from the analysis.
#' @param method the p-value integration method for combining accross multiple treatments. Available methods 
#' are Stouffer, OSP, Fisher, and maxP. The default method is Stouffer.
#' @return a doublet corresponding to the enrichment after integration across all genes or substrates 
#' that belong to the same pathway or kinase, and the size of the mapped genes or substrates to that pathway
#' or kinase.
#' @export
#' 
#' @examples
#' # load the example data
#' data(PM)
#' 
#' # load pathway annotations
#' data(Pathways)
#' 
#' # convert statistics into z-scores
#' PM.zscores <- apply(PM, 2, function(x){qnorm(rank(x)/(nrow(PM)+1))})
#' 
#' # Rotate the matrix by contrast 1, -1, -1 (i.e. up-regulation, down-regulation, dow-regulation).
#' PM.rotated <- rotate3d(PM.zscores, contrast = c(1, -1, -1))
#'
#' # combine rotated statistics across treatments
#' gene.pvalues <- apply(PM.rotated, 1, geneStats)
#' 
#' # compute statistics for all reactome pathways
#' gene.zscores <- qnorm(gene.pvalues, lower.tail = FALSE)
#' gst <- t(sapply(Pathways.reactome, pathwayStats, gene.zscores))
#' 
pathwayStats = function(PGs, T, minSize=5, method="Stouffer"){
  pvalue <- 0
  Z <- T[names(T) %in% PGs]

  if (length(Z) >= minSize) {
    if (method == "Stouffer") {
	   pvalue <- pnorm(sum(Z), 0, sqrt(length(Z)), lower.tail = FALSE)
	} else if (method == "OSP") {
       p <- pnorm(Z, lower.tail = TRUE)
	   pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = TRUE)
	} else if (method == "Fisher") {
	   p <- pnorm(Z, lower.tail = FALSE)
	   pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = FALSE)
    } else if (method == "maxP") {
	   pvalue <- pnorm(max(Z), lower.tail = FALSE)
    }
  } else { 
    pvalue <- NA
  }
  
  result <- list()
  result$pvalue <- pvalue
  result$size <- length(Z)
  return (result)
}
