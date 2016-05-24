#' @title Probability distribution of a CUB model with covariates for both feeling and uncertainty
#' @aliases probcubpq
#' @description Compute the probability distribution of a CUB model with covariates for both the feeling 
#' and the uncertainty components.
#' @export probcubpq
#' @usage probcubpq(m, ordinal, Y, W, bet, gama)
#' @keywords distribution
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param Y Matrix of covariates for explaining the uncertainty component
#' @param W Matrix of covariates for explaining the feeling component
#' @param bet Vector of parameters for the uncertainty component, whose length equals 
#' NCOL(Y) + 1 to include an intercept term in the model (first entry)
#' @param gama Vector of parameters for the feeling component, whose length equals 
#' NCOL(W)+1 to include an intercept term in the model (first entry)
#' @return A vector of the same length as ordinal, whose i-th component is the probability of the 
#' i-th rating according to a CUB distribution with given covariates for both feeling and uncertainty, 
#' and specified coefficients vectors 
#' @seealso \code{\link{bitgama}}, \code{\link{probcub00}}, \code{\link{probcubp0}}, \code{\link{probcub0q}}
#' @references 
#' Piccolo D. (2006). Observed Information Matrix for MUB Models, 
#' \emph{Quaderni di Statistica}, \bold{8}, 33--78 \cr
#' Piccolo D. and D'Elia A. (2008). A new approach for modelling consumers' preferences, \emph{Food Quality and Preference},
#' \bold{18}, 247--259 \cr
#' Iannario M. and Piccolo D. (2012). CUB models: Statistical methods and empirical evidence, in: 
#' Kenett R. S. and Salini S. (eds.), \emph{Modern Analysis of Customer Surveys: with applications using R}, 
#' J. Wiley and Sons, Chichester, 231--258
#' @examples
#' m<-8
#' n<-100
#' ordinal<-sample(1:m,n,replace=TRUE)
#' Y<-sample(c(0,1),n,replace=TRUE)
#' W<-sample(c(0,1,2),n,replace=TRUE)
#' bet<-c(0.3,-0.4)
#' gama<-c(0.1,-1.5)
#' probi<-probcubpq(m, ordinal, Y, W, bet, gama)



probcubpq <-
function(m,ordinal,Y,W,bet,gama){
  logis(Y,bet)*(bitgama(m,ordinal,W,gama)-1/m)+1/m
}
