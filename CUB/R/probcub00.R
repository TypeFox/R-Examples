#' @title Probability distribution of a CUB model without covariates
#' @description Compute the probability distribution of a CUB model without covariates.
#' @aliases probcub00
#' @usage probcub00(m, pai, csi)
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter 
#' @param csi Feeling parameter
#' @export probcub00
#' @return The vector of the probability distribution of a CUB model
#' @keywords distribution
#' @references 
#' Piccolo D. (2003). On the moments of a mixture of uniform and shifted binomial random variables. 
#' \emph{Quaderni di Statistica}, \bold{5}, 85--104\cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#' J. Wiley and Sons, Chichester, 231--258 
#' @examples 
#' m<-9
#' pai<-0.3
#' csi<-0.8
#' pr<-probcub00(m, pai, csi)
#' plot(1:m,pr,type="h",main="CUB probability distribution",xlab="Ordinal categories")
#' points(1:m,pr,pch=19)



probcub00 <-function(m,pai,csi){
  pai*(probbit(m,csi)-1/m)+1/m
  }
