#' @title Normalized Gini heterogeneity index
#' @description Compute the normalized Gini heterogeneity 
#' index for a given discrete probability distribution.
#' @usage gini(m, prob)
#' @aliases gini
#' @param m Number of categories
#' @param prob Vector of the probability distribution or relative frequencies
#' @keywords univar
#' @export gini
#' @seealso \code{\link{laakso}}
#' @examples 
#' m<-7
#' prob<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' gini(m, prob)


gini <-
function(m,prob){
  (1-sum(prob^2))*m/(m-1)}
