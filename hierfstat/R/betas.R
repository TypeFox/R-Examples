#####################################################
#'
#' @title Estimate \eqn{\beta}s per population and a bootstrap confidence interval
#' 
#' @description Estimate \eqn{\beta}s per population and a bootstrap confidence interval 
#' 
#' @usage betas(dat,nboot=0,lim=c(0.025,0.975),diploid=TRUE)
#' 
#' @param dat data frame with genetic data and pop identifier
#' @param nboot number of bootstrap samples
#' @param lim width of the bootstrap confidence interfal
#' @param diploid whether the data comes from a diploid organism
#' 
#' @return betaiovl Average \eqn{\beta_i} over loci
#' @return ci The bootstrap confidence interval
#' @return Hi Within population gene diversities
#' @return Hb Between populations gene diversities
#' 
#' @details  Beware, only valid for non inbred diploids!
#' 
#' @author Jerome Goudet \email{jerome.goudet@@unil.ch}
#'
#'@examples
#' dat<-sim.genot(size=100,N=c(100,1000,10000),nbloc=50,nbal=10)
#' betas(dat,nboot=100)$ci 
#'
#'@export
#'
#####################################################
betas<-function(dat,nboot=0,lim=c(0.025,0.975),diploid=TRUE){
  ratio.Hi.Hb<-function(x){
    dum<-which(!is.na(x))
    sum(x[dum])/sum(Hb[dum]) 
  }
  if (is.genind(dat)) dat<-genind2hierfstat(dat)
  
  pfr<-pop.freq(dat,diploid)
  pfr2<-lapply(pfr,function(x) t(x) %*% x)
  nl<-dim(dat)[2]-1
  np<-length(table(dat[,1]))
  ns<-ind.count.n(dat)
  if (diploid) ns<-ns*2
  ns[is.na(ns)]<-0
  Hi<-matrix(numeric(np*nl),ncol=nl)
  Hb<-numeric(nl)
  for (il in 1:nl){
    Hi[,il]<-ns[,il]/(ns[,il]-1)*(1-diag(pfr2[[il]]))
    #   Hi[is.na(Hi)]<-0.0
    npl<-sum(ns[,il]>0,na.rm=TRUE)     
    Hb[il]<-1-1/npl/(npl-1)*sum((pfr2[[il]]-diag(diag(pfr2[[il]]))),na.rm=TRUE)     
  }
  betai<-1-apply(Hi,1,ratio.Hi.Hb)
  if (nboot<100){
    return(list(betaiovl=betai,Hi=Hi,Hb=Hb))
  }
  else
  {
    if (nl<10) {
      warning("Less than 10 loci, can't estimate Conf. Int.")
      return(list(betaiovl=betai,Hi=Hi,Hb=Hb))
      }
    
    boot.bi<-matrix(numeric(nboot*np),nrow=nboot)
    nls<-apply(ns,1,function(x) which(x>0))
    if (is.matrix(nls)){
    for (ib in 1:nboot){
      for (ip in 1:np){
        dum<-sample(nls[,ip],replace=TRUE)
        boot.bi[ib,ip]<-1-sum(Hi[ip,dum])/sum(Hb[dum])
      }}}
    else{
      for (ib in 1:nboot){
        for (ip in 1:np){
          dum<-sample(nls[[ip]],replace=TRUE)
          boot.bi[ib,ip]<-1-sum(Hi[ip,dum])/sum(Hb[dum])
          
    }}}
    bi.ci<-apply(boot.bi,2,quantile,lim,na.rm=TRUE)
    return(list(betaiovl=betai,ci=bi.ci,Hi=Hi,Hb=Hb))
  }
}
