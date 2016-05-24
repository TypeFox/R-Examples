#' Likelihood of the data and parameters.
#' 
#' Calculates the log-likelihood of the data with the current parameters and
#' Kendall's distance.
#' 
#' 
#' @param z Probability of each cluster membership.
#' @param p Proportion in each cluster.
#' @param C.lam Vector of normalizing coefficients for Mallows' model.
#' @param lambda Current spread parameters
#' @param all.dists.data All distances from the data to the modal sequences.
#' @return Current log-likelihood of the data with the current parameters.
#' @author Erik Gregory
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords likelihood Mallow

Likelihood <-
function(z, p, C.lam, lambda, all.dists.data) {
  to.add <- log(C.lam * p)
  by.lam <- t(-lambda*t(all.dists.data))
  by.lam <- sweep(by.lam, 2, to.add, "+")
  fin <- by.lam * z
  like <- sum(fin)
  return(like)
}
