`HotellingsT.internal` <- function(X,Y=NULL,mu, test)
{
 n<-dim(X)[1]
 p<-dim(X)[2]

 if(is.null(Y))     #one sample case
 {
  test.statistic<-n*as.numeric(t(colMeans(X)-mu)%*%solve(cov(X))%*%(colMeans(X)-mu))*switch(test,f=(n-p)/(p*(n-1)),chi=1)
  df.1<-p
  df.2<-switch(test,f=n-p,chi=NA) 
  p.value<-1-switch(test,f=pf(test.statistic,df.1,df.2),chi=pchisq(test.statistic,df.1))
  return(list(test.statistic=test.statistic,p.value=p.value,df.1=df.1,df.2=df.2))
 }
 
 # else two sample case
 n1<-n
 n2<-dim(Y)[1]
 Xmeans<-colMeans(X)
 Ymeans<-colMeans(Y)
 X.diff<-sweep(X,2,Xmeans)
 Y.diff<-sweep(Y,2,Ymeans)
 S.pooled<-1/(n1+n2-2)*(t(X.diff)%*%X.diff+t(Y.diff)%*%Y.diff)
 test.statistic<-n1*n2/(n1+n2)*t(Xmeans-Ymeans-mu)%*%solve(S.pooled)%*%(Xmeans-Ymeans-mu)*switch(test,f=(n1+n2-p-1)/(p*(n1+n2-2)),chi=1)
 df.1<-p
 df.2<-switch(test,f=n1+n2-p-1,chi=NA)
 p.value<-1-switch(test,f=pf(test.statistic,df.1,df.2),chi=pchisq(test.statistic,df.1))
 list(test.statistic=test.statistic,p.value=p.value,df.1=df.1,df.2=df.2)
}
