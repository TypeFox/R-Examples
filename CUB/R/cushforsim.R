#' @title Simulation routine for CUSH models without covariates
#' @aliases cushforsim
#' @description Fit CUSH models without covariates to given ordinal data. It is useful for simulation 
#' experiments since it performs the same steps as \code{\link{CUSH}}, but with no printed output.
#' @usage cushforsim(m, ordinal, shelter)
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param shelter Category corresponding to the shelter choice
#' @return An object of the class "CUSH", with null output for $BIC since the routine is only for
#'  simulation purposes
#' @seealso \code{\link{CUSH}}
#' @export cushforsim
#' @keywords htest
#' @examples
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,45])
#' shelter<-1
#' simul<-cushforsim(m, ordinal, shelter)
#' simul$estimates
#' simul$loglik
#' simul$varmat



cushforsim <-
function(m,ordinal,shelter){
  freq<-tabulate(ordinal,nbins=m); n<-length(ordinal); aver<-mean(ordinal);
  fc<-freq[shelter]/n
  deltaest<-max(0.01,(m*fc-1)/(m-1))        ### sufficient unbiased estimator of delta
  esdelta<-sqrt((1-deltaest)*(1+(m-1)*deltaest)/(n*(m-1)))
  varmat<-esdelta^2
  wald<-deltaest/esdelta
  loglik<-n*((1-fc)*log(1-deltaest)+fc*log(1+(m-1)*deltaest)-log(m))
  AICCUSH<- -2*loglik+2
  BICCUSH<- -2*loglik+log(n)
  llunif<- -n*log(m); csisb<-(m-aver)/(m-1)
  llsb<-loglikcub00(m,freq,1,csisb)
  nonzero<-which(freq!=0)
  logsat<- -n*log(n)+sum((freq[nonzero])*log(freq[nonzero]))
  devian<-2*(logsat-loglik)
  LRT<-2*(loglik-llunif)
  #####################################################################
  # Assignments as global variables: assign('name',value,pos=1)
  #####################################################################
  #   assign('deltaest',deltaest,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('loglik',loglik,pos=1)
  results<-list('estimates'=round(deltaest,digits=5), 'loglik'=loglik,'varmat'=varmat)
  
}
