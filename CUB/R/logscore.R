#' @title Logarithmic score
#' @description Compute the logarithmic score of a CUB model with covariates both for the uncertainty 
#' and the feeling parameters.
#' @aliases logscore
#' @export logscore
#' @usage logscore(m, ordinal, Y, W, bet, gama)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of covariates for explaining the uncertainty component
#' @param W Matrix of covariates for explaining the feeling component
#' @param bet Vector of parameters for the uncertainty component, with length NCOL(Y)+1 
#' to account for an intercept term (first entry of bet)
#' @param gama Vector of parameters for the feeling component, with length NCOL(W)+1 
#' to account for an intercept term (first entry of gama)
#' @references 
#' Tutz, G. (2012). \emph{Regression for Categorical Data}, Cambridge University Press, Cambridge
#' @keywords htest
#' @examples
#' data(relgoods)
#' m<-10
#' gender<-relgoods[,2]
#' smoking<-relgoods[,12]
#' ordinal<-relgoods[,40]
#' nona<-na.omit(cbind(ordinal,gender,smoking))
#' ordinal<-nona[,1]
#' gender<-nona[,2]
#' smoking<-nona[,3]
#' bet<-c(-0.45,-0.48)
#' gama<-c(-0.55,-0.43)
#' logscore(m, ordinal, Y=smoking, W=gender, bet, gama)


logscore <-
function(m,ordinal,Y,W,bet,gama){
  pr<-probcubpq(m,ordinal,Y,W,bet,gama)
  return(-2*sum(log(pr)))
}
