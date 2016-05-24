#' Update Proportion in each cluster.
#' 
#' Updates the proportion of data assigned to each cluster.
#' 
#' 
#' @param z Probabilities that each sequence is in each cluster.
#' @return Proportion of data in each cluster.
#' @author Erik Gregory
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#' @keywords proportion
UpdateP <-
function(z) {
  p <- colSums(z)
  p <- p/nrow(z)
  return(p)
}
