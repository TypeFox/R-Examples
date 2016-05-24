#' @title Log-likelihood function of CUBE models for ordinal data 
#' @aliases loglikcuben
#' @description Compute the log-likelihood function of a CUBE model 
#' without covariates for ordinal responses, possibly with different
#' vectors of parameters for each observation.
#' @usage loglikcuben(m, ordinal, assepai, assecsi, assephi)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param assepai Vector of uncertainty parameters for the given 
#' observations (with the same length as ordinal)
#' @param assecsi Vector of feeling parameters for the given observations
#' (with the same length as ordinal)
#' @param assephi Vector of overdispersion parameters for the given 
#' observations (with the same length as ordinal)
#' @export loglikcuben
#' @seealso \code{\link{loglikCUBE}}
#' @keywords htest
#' @examples
#' m<-8
#' n0<-230;  n1<-270;
#' bet<-c(-1.5,1.2)
#' gama<-c(0.5,-1.2)
#' alpha<-c(-1.2,-0.5)
#' pai0<-1/(1+exp(-bet[1])); csi0<-1/(1+exp(-gama[1])); phi0<-exp(alpha[1]);
#' ordinal0<-simcube(n0,m,pai0,csi0,phi0)
#' pai1<-1/(1+exp(-sum(bet))); csi1<-1/(1+exp(-sum(gama))); phi1<-exp(sum(alpha));
#' ordinal1<-simcube(n1,m,pai1,csi1,phi1)
#' ordinal<-c(ordinal0,ordinal1)
#' assepai<-c(rep(pai0,n0),rep(pai1,n1))
#' assecsi<-c(rep(csi0,n0),rep(csi1,n1))
#' assephi<-c(rep(phi0,n0),rep(phi1,n1))
#' lli<-loglikcuben(m,ordinal,assepai,assecsi,assephi)


loglikcuben <-
function(m,ordinal,assepai,assecsi,assephi){
  n<-length(ordinal)
  prob<-matrix(NA,nrow=n,ncol=m)
  pconi<-rep(NA,n)
  for (i in 1:n){
    prob[i,]<-t(probcube(m,assepai[i],assecsi[i],assephi[i]))
    pconi[i]<-prob[i,ordinal[i]]
  }
   return(sum(log(pconi)))
}
