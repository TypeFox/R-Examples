#' @export
resVar<-function(x, adj=c("OLS", "ML")){
  if(!inherits(x, "setar"))
    stop("Object must be created by setar()")
  adj<-match.arg(adj)
  obj<-x$model.specific
  z<-obj$z
  regime<-tail(obj$regime, nrow(x$str$xx))
  sigGen<-crossprod(na.omit(x$res))/(length(x$res)-ifelse(adj=="OLS",obj$k,0))
  names(sigGen)<-"Total"
  nc<-ifelse(adj=="OLS", length(obj$IncNames), 0)
  ML<-ifelse(adj=="OLS", length(obj$ML), 0)
  MH<-ifelse(adj=="OLS", length(obj$MH), 0)
  MM<-ifelse(adj=="OLS", length(obj$MM), 0)

  if(obj$nthresh==1){
    isL <- ifelse(regime==ifelse(obj$restriction=="none",1,2), 1,0)
    isH <- 1-isL
    resL<-isL*x$res
    sigL<-crossprod(na.omit(resL))/(sum(isL)-ML-nc)
    resH<-isH*x$res
    sigH<-crossprod(na.omit(resH))/(sum(isH)-MH-nc)
    sigReg<-c(sigL, sigH)
    names(sigReg)<-c("L", "H")
  }
  if(obj$nthresh==2){
    isL<-ifelse(regime==1, 1,0)
    isM<-ifelse(regime==2, 1,0)
    isH<-1-isL-isM
    resL<-isL*x$res
    resM<-isM*x$res
    resH<-isH*x$res
    sigL<-crossprod(na.omit(resL))/(sum(isL)-ML-nc)
    sigM<-crossprod(na.omit(resM))/(sum(isM)-MM-nc)
    sigH<-crossprod(na.omit(resH))/(sum(isH)-MH-nc)
    sigReg<-c(sigL, sigM, sigH)
    names(sigReg)<-c("L", "M","H")
  }
c(sigGen,sigReg)
}

## Method for nlVar, in case resVar made once as method
## Problems now: use n-1 instead of either n or n-k
resVar.nlVar <- function(x){
  
  nobs <- x$nobs*x$t  
  
  ##extract regime and residuals:
  resid <- residuals(x)
  regs <- regime(x, initVal=FALSE, timeAttr=FALSE)
  
  vars <- aggregate(x=resid, list(regs), var)
  colnames(vars)[1]<- "Regime"
  
  vars
}


if(FALSE){
library(tsDyn)
environment(resVar)<-environment(star)
a<-setar(lynx, m=2)

resVar(a)

aa<-setar(lynx, m=2, nthresh=2)
aa
resVar(aa)

}

