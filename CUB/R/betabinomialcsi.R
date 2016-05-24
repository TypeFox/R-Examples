#' @title Beta-Binomial probabilities of ordinal responses, given feeling parameter for each observation
#' @description Compute the Beta-Binomial probabilities of the given ordinal responses, with feeling 
#' parameter specified for each observation, 
#' and with the same overdispersion parameter for all the responses.
#' @aliases betabinomialcsi
#' @usage betabinomialcsi(m, ordinal, csivett, phi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param csivett  Vector of feeling parameters of the Beta-Binomial distribution for the given ordinal 
#' responses
#' @param phi Overdispersion parameter of the Beta-Binomial distribution 
#' @export betabinomialcsi
#' @return A vector of the same length as ordinal, containing the Beta-Binomial probability of each observation 
#' for the corresponding feeling parameter and for the specified overdispersion parameter
#' @seealso   \code{\link{betar}}, \code{\link{betabinomial}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data,
#'  \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
#' @keywords distribution
#' @examples 
#' data(relgoods)
#' m<-10
#' ordinal<-relgoods[,37]
#' age<-2014-relgoods[,4]
#' lage<-log(age)-mean(log(age))
#' nona<-na.omit(cbind(ordinal,lage))
#' ordinal<-nona[,1]
#' W<-nona[,2]
#' gama<-c(-0.61,-0.31)
#' phi<-0.16 
#' csivett<-logis(W,gama)
#' pr<-betabinomialcsi(m,ordinal,csivett,phi)


betabinomialcsi <-function(m,ordinal,csivett,phi){
  n<-length(ordinal)
  betabin<-rep(NA,m)
  for(i in 1:n){
    bebeta<-betar(m,csivett[i],phi)
    betabin[i]<-bebeta[ordinal[i]]   
  }
  return(betabin)
}


