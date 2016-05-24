"mkle.ci" <- function(data, bw=2*sd(data), alpha=0.1, kernel=c("gaussian", "epanechnikov", "rectangular", 
                      "triangular", "biweight", "cosine", "optcosine"), method=c("percentile", 
                      "wald","boott"), B=1000,gridsize=2^14){
 
  method <- match.arg(method,several.ok=TRUE)
  if(alpha>=0.5 || alpha <=0 )stop('alpha must lie between 0 and 0.5')
  out<-NULL
  if (any(method=='boott')) {
    n<-length(data)
    bootdata<-matrix(sample(data,length(data)*B,replace=TRUE),ncol=B)
    jackmkles<-matrix(NA,ncol=n,nrow=B)
    for(i in 0:(n-1)){
      jackdata<-apply(bootdata,2,subset,c(rep(TRUE,i),FALSE,rep(TRUE,n-i-1)))[1:(n-1),]
      jackmkles[,i+1]<-apply(jackdata,2,mkle,bw=bw,gridsize=gridsize,kernel=kernel)
    }
    sestars<-sqrt(apply(jackmkles,1,var)*(n-1)^2/n)
    mkles<-apply(bootdata,2,mkle,bw=bw,gridsize=gridsize,kernel=kernel)
    midpoint<-mkle(data,bw=bw,gridsize=gridsize,kernel=kernel)
    tstar<-quantile((mkles - midpoint) / sestars, c(1-alpha/2,alpha/2), type=5)
    out<-rbind(out,data.frame(method='boott',lower=midpoint-tstar[1]*sd(mkles),upper=midpoint-tstar[2]*sd(mkles)))
  }else{
    mkles<-apply(matrix(sample(data,length(data)*B,replace=TRUE),ncol=B),2,mkle,bw=bw,gridsize=gridsize,kernel=kernel)
    midpoint<-mkle(data,bw=bw,gridsize=gridsize,kernel=kernel)
  }
  if (any(method=='percentile')){
    out<-rbind(out,data.frame(method='percentile',lower=quantile(mkles,alpha/2,type=5),upper=quantile(mkles,1-alpha/2,type=5)))
  }
  if (any(method=='wald')){
    out<-rbind(out,data.frame(method='wald',lower=midpoint+qnorm(alpha/2)*sd(mkles),upper=midpoint+qnorm(1-alpha/2)*sd(mkles)))
  }

  rownames(out)<-1:length(out[,1])
  return(out)
}
