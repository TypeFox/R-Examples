#' @title Probability distribution of a CUB model with covariates for the feeling component
#' @aliases probcub0q
#' @description Compute the probability distribution of a CUB model with covariates
#'  for the feeling component.
#' @export probcub0q
#' @usage probcub0q(m, ordinal, W, pai, gama)
#' @keywords distribution
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param W Matrix of covariates for explaining the feeling component
#' NCOL(Y)+1 to include an intercept term in the model (first entry)
#' @param pai Uncertainty parameter
#' @param gama Vector of parameters for the feeling component, whose length equals 
#' NCOL(W)+1 to include an intercept term in the model (first entry)
#' @return A vector of the same length as ordinal, whose i-th component is the
#' probability of the i-th observation according to a CUB model with the corresponding values 
#' of the covariates for the feeling component, and coefficients specified in gama
#' @seealso \code{\link{bitgama}}, \code{\link{probcub00}}, \code{\link{probcubp0}}, 
#' \code{\link{probcubpq}}
#' @references 
#' Piccolo D. (2006). Observed Information Matrix for MUB Models, 
#' \emph{Quaderni di Statistica}, \bold{8}, 33--78 \cr
#' Piccolo D. and D'Elia A. (2008). A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
#' \bold{18}, 247--259 \cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#' J. Wiley and Sons, Chichester, 231--258
#' @examples
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,29]
#' gender<-relgoods[,2]
#' data<-na.omit(cbind(ordinal,gender))
#' ordinalnew<-data[,1]
#' W<-data[,2]
#' pai<-0.44
#' gama<-c(-0.91,-0.7)
#' pr<-probcub0q(m,ordinalnew,W,pai,gama)

probcub0q <-
function(m,ordinal,W,pai,gama){
  pai*(bitgama(m,ordinal,W,gama)-1/m)+1/m
}
