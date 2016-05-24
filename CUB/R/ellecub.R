#' @title Log-likelihood function of a CUB model without covariates
#' @description Compute the log-likelihood function of a CUB model without covariates fitting 
#' ordinal responses, possibly with different vector of parameters for each observation.
#' @aliases ellecub
#' @usage ellecub(m, ordinal, assepai, assecsi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param assepai Vector of uncertainty parameters for given observations 
#' (with the same length as ordinal)
#' @param assecsi Vector of feeling parameters for given observations
#'  (with the same length as ordinal)
#' @seealso \code{\link{loglikCUB}}
#' @keywords htest
#' @export ellecub
#' @examples
#' m<-7
#' n0<-230
#' n1<-270
#' bet<-c(-1.5,1.2)
#' gama<-c(0.5,-1.2)
#' pai0<-1/(1+exp(-bet[1])); csi0<-1/(1+exp(-gama[1]));
#' pai1<-1/(1+exp(-sum(bet))); csi1<-1/(1+exp(-sum(gama)));
#' ordinal0<-simcub(n0,m,pai0,csi0)
#' ordinal1<-simcub(n1,m,pai1,csi1)
#' ordinal<-c(ordinal0,ordinal1)
#' assepai<-c(rep(pai0,n0),rep(pai1,n1))
#' assecsi<-c(rep(csi0,n0),rep(csi1,n1))
#' lli<-ellecub(m,ordinal,assepai,assecsi)



ellecub <-
function(m,ordinal,assepai,assecsi){
  prob<-assepai*(dbinom(0:(m-1),m-1,1-assecsi)-1/m)+1/m
  pconi<-prob[ordinal]
  return(sum(log(pconi)))
}
