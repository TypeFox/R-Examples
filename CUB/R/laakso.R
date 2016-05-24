#' @title Normalized Laakso and Taagepera heterogeneity index
#' @description Compute the normalized Laakso & Taagepera heterogeneity index for a given 
#' discrete probability distribution.
#' @aliases laakso
#' @export laakso
#' @usage laakso(m, prob)
#' @param m Number of categories
#' @param prob Vector of length \eqn{m} of a probability or relative frequency distribution
#' @seealso \code{\link{gini}}
#' @keywords univar
#' @examples
#' m<-7
#' prob<-c(0.04,0.04,0.05,0.10,0.21,0.32,0.24)
#' laakso(m,prob)



laakso <-
function(m,prob){
  1/(m/gini(m,prob)-m+1)}
