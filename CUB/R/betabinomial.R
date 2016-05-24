#' @title Beta-Binomial probabilities of ordinal responses, with feeling and overdispersion parameters
#' for each observation
#' @description Compute the Beta-Binomial probabilities of ordinal responses, given feeling and overdispersion
#' parameters for each observation.
#' @aliases betabinomial
#' @usage betabinomial(m, ordinal, csivett, phivett)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param csivett  Vector of feeling parameters of the Beta-Binomial distribution for the given ordinal responses
#' @param phivett Vector of overdispersion parameters of the Beta-Binomial distribution for the given ordinal 
#' responses
#' @export betabinomial
#' @return A vector of the same length as ordinal, containing the Beta-Binomial probability of each observation,
#'  for the corresponding feeling and overdispersion parameters.
#' @details The Beta-Binomial distribution is the Binomial distribution in which the probability of success at
#'  each trial is not fixed but random and follows the Beta distribution. It is frequently used in Bayesian 
#'  statistics, empirical Bayes methods and classical statistics as an overdispersed binomial distribution. 
#' @seealso  \code{\link{betar}}, \code{\link{betabinomialcsi}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
#' @keywords distribution
#' @examples 
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,37]
#' age<-2014-relgoods[,4]
#' lage<-log(age)-mean(log(age))
#' nona<-na.omit(cbind(ordinal,lage))
#' ordinal<-nona[,1]
#' gama<-c(-0.6, -0.3)
#' csivett<-logis(lage,gama)
#' alpha<-c(-2.3,0.92)
#' phivett<-1/(-1+ 1/(logis(lage,alpha)))
#' pr<-betabinomial(m, ordinal, csivett, phivett)



betabinomial <-
function(m,ordinal,csivett,phivett){
  n<-length(ordinal)
  betabin<-rep(NA,n)
  for(i in 1:n){
    bebeta<-betar(m,csivett[i],phivett[i])
    betabin[i]<-bebeta[ordinal[i]]   
  }
  return(betabin)
}

