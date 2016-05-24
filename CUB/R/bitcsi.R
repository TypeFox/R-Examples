#' @title Shifted Binomial probabilities of ordinal responses
#' @description Compute the shifted Binomial probabilities of ordinal responses.
#' @aliases bitcsi
#' @usage bitcsi(m, ordinal, csi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param csi  Feeling parameter of the shifted Binomial distribution
#' @export bitcsi
#' @return A vector of the same length as ordinal, where each entry is the shifted Binomial probability 
#' of the corresponding observation 
#' @seealso  \code{\link{probcub00}}, \code{\link{probcubp0}},  \code{\link{probcub0q}}
#' @references Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables,
#'  \emph{Quaderni di Statistica}, \bold{5}, 85--104
#' @keywords distribution
#' @examples 
#' data(univer)
#' m<-7
#' ordinal<-univer[,8]
#' csi<-0.7
#' pr<-bitcsi(m,ordinal,csi)


bitcsi <-function(m,ordinal,csi){
  base<-log(1-csi)-log(csi)
  const<-exp(m*log(csi)-log(1-csi))
  const*kkk(m,ordinal)*exp(base*ordinal)
}

