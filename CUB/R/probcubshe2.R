#' @title probcubshe2
#' @aliases probcubshe2
#' @description Probability distribution of a CUB model with explicit shelter effect
#' @usage probcubshe2(m, pai, csi, delta, shelter)
#' @export probcubshe2
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @param delta Shelter parameter
#' @param shelter Category corresponding to the shelter choice
#' @return The vector of the probability distribution of a CUB model with explicit shelter effect
#' @details A CUB model with explicit shelter effect is a mixture of two components: 
#' a CUB distribution with uncertainty parameter \eqn{\pi}  and feeling parameter \eqn{\xi}, 
#' and  a degenerate distribution with unit mass at the shelter category ("shelter"), 
#' with mixing coefficient specified by \eqn{\delta}.
#' @references 
#' Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
#' \emph{Statistical Methods and Applications}, \bold{21}, 1--22 \cr
#' Iannario M. and Piccolo D. (2014). Inference for CUB models: a program in R,
#' \emph{Statistica & Applicazioni}, \bold{XII} n.2, 177--204 \cr
#' Iannario M. and Piccolo D. (2014). A comprehensive approach to ordinal data modelling,
#' \emph{Working paper}
#' @seealso \code{\link{probcubshe1}}, \code{\link{probcubshe3}}
#' @keywords distribution 
#' @examples
#' m<-8
#' pai1<-0.5
#' pai2<-0.3
#' csi<-0.4
#' shelter<-6
#' delta<-1-pai1-pai2
#' pai<-pai1/(1-delta)
#' pr2<-probcubshe2(m, pai, csi, delta, shelter)
#' plot(1:m,pr2,type="h", main="CUB probability distribution with 
#' explicit shelter effect",xlab="Ordinal categories")
#' points(1:m,pr2,pch=19)

probcubshe2 <-
function(m,pai,csi,delta,shelter)
{delta*ifelse(seq(1,m)==shelter,1,0)+(1-delta)*(pai*probbit(m,csi)+(1-pai)*(1/m))}
