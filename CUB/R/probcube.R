#' @title Probability distribution of a CUBE model without covariates
#' @aliases probcube
#' @description Compute the probability distribution of a CUBE model without covariates.
#' @usage probcube(m, pai, csi, phi)
#' @export probcube
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @param phi Overdispersion parameter
#' @return The vector of the probability distribution of a CUBE model without covariates
#' @seealso \code{\link{betar}}, \code{\link{betabinomial}}
#' @references 
#' Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786 
#' @keywords distribution
#' @examples 
#' m<-9
#' pai<-0.3
#' csi<-0.8
#' phi<-0.1
#' pr<-probcube(m, pai, csi, phi)
#' plot(1:m,pr,type="h", main="CUBE probability distribution",xlab="Ordinal categories")
#' points(1:m,pr,pch=19)

probcube <-
function(m,pai,csi,phi){
  pai*(betar(m,csi,phi)-1/m)+1/m
}
