#' @title probcubshe3
#' @aliases probcubshe3
#' @description Probability distribution of a CUB model with explicit shelter effect: 
#' satisficing interpretation
#' @usage probcubshe3(m, lambda, eta, csi, shelter)
#' @export probcubshe3
#' @keywords distribution
#' @param m Number of ordinal categories
#' @param lambda Mixing coefficient for the shifted Binomial component
#' @param eta Mixing coefficient for the mixture of the uncertainty component and the 
#' shelter effect
#' @param csi Feeling parameter
#' @param shelter Category corresponding to the shelter choice
#' @return The vector of the probability distribution of a CUB model with shelter effect
#' @details The "satisficing interpretation" provides a parametrization for CUB models with explicit 
#' shelter effect as a mixture of two components: a shifted Binomial distribution with feeling parameter
#' \eqn{\xi} (meditated choice), and a mixture of a degenerate distribution with unit mass at the shelter
#' category ("shelter") and a discrete uniform distribution over \eqn{m} categories, with mixing 
#' coefficient specified by \eqn{\eta} (lazy selection of a category). Both components of the mixtures
#' are weighted by \eqn{\lambda} coefficient.
#' @references 
#' Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
#' \emph{Statistical Methods and Applications}, \bold{21}, 1--22 \cr
#' Iannario M. and Piccolo D. (2014). Inference for CUB models: a program in R,
#' \emph{Statistica & Applicazioni}, \bold{XII} n.2, 177--204 \cr
#' Iannario M. and Piccolo D. (2014). A comprehensive approach to ordinal data modelling,
#' \emph{Working paper}
#' @seealso \code{\link{probcubshe1}}, \code{\link{probcubshe2}}
#' @examples
#' m<-8
#' pai1<-0.5
#' pai2<-0.3
#' csi<-0.4
#' shelter<-6
#' lambda<-pai1
#' eta<-1-pai2/(1-pai1)
#' pr3<-probcubshe3(m, lambda, eta, csi, shelter)
#' plot(1:m,pr3,type="h", main="CUB probability distribution with explicit 
#' shelter effect",xlab="Ordinal categories")
#' points(1:m,pr3,pch=19)



probcubshe3 <-
function(m,lambda,eta,csi,shelter){
  lambda*probbit(m,csi)+(1-lambda)*((1-eta)/m  + eta*ifelse(seq(1,m)==shelter,1,0))
}
