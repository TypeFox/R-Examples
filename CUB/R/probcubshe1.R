#' @title probcubshe1
#' @aliases probcubshe1
#' @description Probability distribution of an extended CUB model with a shelter effect.
#' @usage probcubshe1(m, pai1, pai2, csi, shelter)
#' @export probcubshe1
#' @keywords distribution
#' @param m Number of ordinal categories
#' @param pai1 Mixing coefficient for the shifted Binomial component of the mixture distribution
#' @param pai2 Mixing coefficient for the discrete Uniform component of the mixture distribution
#' @param csi Feeling parameter
#' @param shelter Category corresponding to the shelter choice
#' @return The vector of the probability distribution of an extended CUB model with a shelter effect 
#' at the shelter category
#' @details An extended CUB model is a mixture of three components: a shifted Binomial distribution 
#' with probability of success \eqn{\xi}, a discrete uniform distribution with support \eqn{\{1,...,m\}},
#'  and a degenerate distribution with unit mass at the shelter category ("shelter").
#' @references 
#' Iannario M. (2012). Modelling \emph{shelter} choices in a class of mixture models for ordinal responses,  
#' \emph{Statistical Methods and Applications}, \bold{21}, 1--22 \cr
#' Iannario M. and Piccolo D. (2014). Inference for CUB models: a program in R,
#'  \emph{Statistica & Applicazioni}, \bold{XII} n.2, 177--204 \cr
#'  Iannario M. and Piccolo D. (2014). A comprehensive approach to ordinal data modelling,
#'   \emph{Working paper}
#' @seealso \code{\link{probcubshe2}}, \code{\link{probcubshe3}}
#' @examples
#' m<-8
#' pai1<-0.5
#' pai2<-0.3
#' csi<-0.4
#' shelter<-6
#' pr<-probcubshe1(m, pai1, pai2, csi, shelter)
#' plot(1:m,pr,type="h", main="Extended CUB probability distribution with shelter effect",
#' xlab="Ordinal categories")
#' points(1:m,pr,pch=19)

probcubshe1 <-
function(m,pai1,pai2,csi,shelter)
{pai1*probbit(m,csi)+pai2*(1/m)+(1-pai1-pai2)*ifelse(seq(1,m)==shelter,1,0)}
