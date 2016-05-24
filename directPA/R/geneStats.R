#' Molecule Level Statistics
#' 
#' Takes a vector of statistics with each element corresponds to a treatment vs control comparison, 
#' and calculates a combined statistics accross multiple treatments.
#' 
#' @usage geneStats(T, method="OSP")
#' @param T a vector of statistics (z-scores converted) with each element correspond to a treatment vs control 
#' comparison.
#' @param method the p-value integration method for combining accross multiple treatments. Available methods 
#' are Stouffer, OSP, Fisher, and maxP. The default method is OSP.
#' @return a p-value after integration across treatments.
#' @export
#' 
#' @examples
#' # load the example data
#' data(PM)
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
geneStats <- function(T, method="OSP") {
   pvalue <- 0

   if (method == "Stouffer") {
      pvalue <- pnorm(sum(T), 0, sqrt(length(T)), lower.tail =FALSE)
   } else if (method == "OSP") {
      p <- pnorm(T, lower.tail = TRUE)
	    pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = TRUE)
   } else if (method == "Fisher") {
	    p <- pnorm(T, lower.tail = FALSE)
	    pvalue <- pchisq(-2*sum(log(p)), 2*length(p), lower.tail = FALSE)
   } else if (method == "maxP") {
	    pvalue <- pnorm(max(T), lower.tail = FALSE)
   }

   return (pvalue)
}
