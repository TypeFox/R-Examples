#' @title Probability distribution of a CUB model with covariates for the uncertainty component
#' @aliases probcubp0
#' @description Compute the probability distribution of a CUB model with covariates for the 
#' uncertainty component.
#' @export probcubp0
#' @usage probcubp0(m, ordinal, Y, bet, csi)
#' @keywords distribution
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of covariates for explaining the uncertainty component
#' @param bet Vector of parameters for the uncertainty component, whose length equals 
#'  NCOL(Y) + 1 to include an intercept term in the model (first entry)
#' @param csi Feeling parameter
#' @return A vector of the same length as ordinal, whose i-th component is the probability of the i-th 
#' observation according to a CUB model with the corresponding values of the covariates for the 
#' uncertainty component, and coefficients for the covariates specified in bet 
#' @references 
#' Piccolo D. (2006). Observed Information Matrix for MUB Models, 
#' \emph{Quaderni di Statistica}, \bold{8}, 33--78 \cr
#' Piccolo D. and D'Elia A. (2008). A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
#' \bold{18}, 247--259 \cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#' J. Wiley and Sons, Chichester, 231--258
#' @seealso \code{\link{bitgama}}, \code{\link{probcub00}}, \code{\link{probcubpq}}, \code{\link{probcub0q}}
#' @examples
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,29]
#' gender<-relgoods[,2]
#' nona<-na.omit(cbind(ordinal,gender))
#' ordinalnew<-nona[,1]
#' Y<-nona[,2]
#' bet<-c(-0.81,  0.93)
#' csi<-0.20
#' probi<-probcubp0(m,ordinalnew,Y,bet,csi)

probcubp0 <-
function(m,ordinal,Y,bet,csi){
  logis(Y,bet)*(bitcsi(m,ordinal,csi)-1/m)+1/m
}
