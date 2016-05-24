GenomicControl<-function(x, snp.sel)
 {
  if(!inherits(x,"WGassociation"))
   stop("x must be an object of class 'WGassociation'")

  WGchisq<-function(x, model) {
     df<-ifelse(model=="codominant", 2, 1)
     qchisq(x,df, lower.tail=FALSE)
  } 

  if (missing(snp.sel)) snp.sel<-rep(TRUE,nrow(x))
  p<-pvalues(x)[snp.sel,-1]
  chisq.obs<-sapply(1:ncol(p) ,function(x) WGchisq(p[,x],names(p)[x]))
   
  lambda<-apply(chisq.obs,2,median,na.rm=TRUE) 
  names(lambda)<-names(x)[-1]
  
  den<-rep(0.456,ncol(x)-1)
  den[names(lambda)=="codominant"]<-1.388
  lambda<- lambda/den
  
  lambdaOK<-ifelse(lambda<1,1,lambda)
  chisq.corrrected<-sweep(chisq.obs, 2, lambdaOK,FUN="/")
  pOK<-1-pchisq(chisq.corrrected,1)
  pOK[pOK==0]<-NA

  k<-length(names(x))
  attr(x,"pvalues")[,2:k]<-pOK

  cat("\nlambda:\n")
  print(lambda)
  x

}
