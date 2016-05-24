
#' @keywords internal


varcovcubshe <-function(m,pai1,pai2,csi,shelter,n){
  pr<-probcubshe1(m,pai1,pai2,csi,shelter)
  dd<-rep(0,m);dd[shelter]<-1;
  bb<-probbit(m,csi)
  ########################
  aaa<-bb-dd
  bbb<-(1/m)-dd
  c4<-pai1*bb*(m-(1:m)-csi*(m-1))/(csi*(1-csi))
  atilde<-aaa/pr;  btilde<-bbb/pr;  ctilde<-c4/pr;
  
  d11<-sum(aaa*atilde);  d22<-sum(bbb*btilde);   dxx<-sum(c4*ctilde);
  d12<-sum(bbb*atilde);  d1x<-sum(c4*atilde);    d2x<-sum(c4*btilde);
  
  ### Information matrix 
  matinf<-matrix(c(d11,d12,d1x,d12,d22,d2x,d1x,d2x,dxx),nrow=3,byrow=T) 
  ### Var-covar matrix 
  
  
  if(any(is.na(matinf))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=3,ncol=3)
  } else {
    if(det(matinf)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=3,ncol=3)
    } else {
      varmat<-solve(matinf)/n
    }
  }
  
  
  return(varmat)
}
