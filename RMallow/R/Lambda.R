#' Objective function to determine lambda.
#' 
#' Objective function to find the root of in calculating the lambda parameters
#' for each cluster.
#' 
#' 
#' @param lambda lambda value to calculate the function output at.
#' @param rhs Right-hand side of the equation in the referenced paper.
#' @param dists Not used.
#' @param dists.table Table of distances between each sequence and the modal
#' sequence in N! space.
#' @return Output of the objective function to determine the root of.  Goal is
#' zero.
#' @author Erik Gregory
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords Mallow lambda

Lambda <-
function(lambda, rhs, dists, dists.table = NULL) {
  C <- C_lam(lambda, dists.table = dists.table)
  # Can we optimize this using previously-known lambda values?
  out <- sum(as.numeric(names(dists.table))*dists.table*exp(-lambda*as.numeric(names(dists.table)))) - rhs
  return(out)
}
