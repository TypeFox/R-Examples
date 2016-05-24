#' @title Simulation routine for CUB models without covariates
#' @aliases cubforsim
#' @description Fit CUB models without covariates to given ordinal data. It is useful for simulation experiments since it performs the same
#' steps as \code{\link{CUB}} function of the package, but with no printed output.
#' @usage cubforsim(m, ordinal, maxiter = 500, toler = 1e-06)
#' @export cubforsim
#' @keywords htest
#' @param m Number of ordinal categories
#' @param ordinal Vector of ordinal responses
#' @param maxiter Maximum number of iterations allowed for running the optimization algorithm (defaul: maxiter = 500)
#' @param toler Fixed error tolerance for final estimates (default: toler = 1e-6)
#' @return An object of the class "CUB", with null output for $BIC since the routine is only for simulation purposes
#' @import stats
#' @seealso \code{\link{CUB}}, \code{\link{loglikCUB}}
#' @examples
#' data(relgoods)
#' m<-10
#' ordinal<-na.omit(relgoods[,37])
#' simul<-cubforsim(m,ordinal,maxiter=500,toler=1e-6)
#' simul$estimates      # Estimated parameters vector (pai,csi)
#' ###############
#' data(univer)
#' m<-7
#' ordinal<-univer[,12]
#' simul<-cubforsim(m,ordinal)
#' param<-simul$estimates   # Estimated parameters vector (pai,csi)
#' ###############
#' m<-9; n<-500;
#' pai<-0.7
#' csi<-0.4
#' ordinal<-simcub(n,m,pai,csi)
#' simul<-cubforsim(m,ordinal)
#' param<-simul$estimates
#' maxlik<-simul$loglik
#' niter<-simul$niter
#' varmat<-simul$varmat 


cubforsim <-
function(m,ordinal,maxiter=500,toler=1e-6){
  serie<-1:m; freq<-tabulate(ordinal,nbins=m); n<-sum(freq);  
  inipaicsi<-inibest(m,freq); pai<-inipaicsi[1]; csi<-inipaicsi[2];
  loglik<-loglikcub00(m,freq,pai,csi)
  # ********************************************************************
  # ************* E-M algorithm for CUB(0,0) ***************************
  # ********************************************************************
  nniter<-1
  while(nniter<=maxiter){
    likold<-loglik
    bb<-probbit(m,csi)
    aa<-(1-pai)/(m*pai*bb)
    tau<-1/(1+aa)
    ft<-freq*tau
    averpo<-(t(serie)%*%ft)/sum(ft)
    pai<-(t(freq)%*%tau)/n  # updated pai estimate
    csi<-(m-averpo)/(m-1)   # updated csi estimate
    if(csi<0.001){
      csi<-0.01; nniter<-maxiter-1;
    }
    # print(c(pai,csi));
    loglik<-loglikcub00(m,freq,pai,csi)
    liknew<-loglik
    testll<-abs(liknew-likold) ###### print(testll); 
    # OPTIONAL printing: print(cbind(nniter,testll,pai,csi));
    if(testll<=toler) break else {loglik<-liknew}
    # OPTIONAL printing: print(loglik);
    nniter<-nniter+1
  }
  ######
  if(csi>0.999) csi<-0.99         ### to avoid division by 0 !!!
  if(csi<0.001) csi<-0.01         ### to avoid division by 0 !!!
  if(pai<0.001) pai<-0.01         ### to ensure identifiability !!!
  ################################################################
  varmat<-varcovcub00(m,ordinal,pai,csi);
  theorpr<-probcub00(m,pai,csi)
  dissimi<-dissim(theorpr,freq/n)
  ################################################################
  #   assign('loglik',loglik,pos=1)
  #   assign('pai',pai,pos=1)
  #   assign('csi',csi,pos=1)
  #   assign('varmat',varmat,pos=1)
  #   assign('nniter',nniter,pos=1)
  stime<-c(pai,csi)
  results<-list('estimates'=round(stime,digits=5), 'loglik'=loglik, 'niter'=nniter,'varmat'=varmat)
  
}
