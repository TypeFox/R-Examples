#' Calculate the Weights for Bubble Plots
#' 
#' Calculate weights for the bubbles in the plot method of \code{cds} objects. The relative 
#' frequencies within a dset of groups are used to calculate the size of the bubble so that the area
#' of the bubble is proportional to the relative frequency of the rating category within that group.
#' 
#' @param dat A data set from which to derive the relative frequencies
#' @param grp A vector giving the group memberships.
#' @param q An integer such that the rating scale is \code{1:q}.
#' @param fact A schrinkage factor.
#' @author Pieter Schoonees
#' @keywords utility
#' @export 
calc.wt.bubbles <- function(dat, grp, q, fact = 0.12) {
  dat.split <- split(dat, grp)
  K <- length(unique(grp))
  counts <- t(sapply(dat.split, function(x) table(factor(unlist(x), levels = 1:q))))
  out <- fact * sqrt(counts/sum(counts))
#   out <- fact*t(counts/matrix(colSums(counts), nrow = q, 
#                               ncol = K, byrow = TRUE))
#   out <- fact*t(counts/matrix(apply(counts, 2, use.fun[[1]]), nrow = q, 
#                               ncol = K, byrow = TRUE))
  out
}
