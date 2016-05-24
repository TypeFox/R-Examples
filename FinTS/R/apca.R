"apca" <- function(x, nf){
# This is a simple program to perform
# Asymptotic Principal Component Analysis.
# The number of nf must be given. 
# The selection criteria will be added later.
#  Written by Ruey S. Tsay
#
  x1=as.matrix(x)
  T.=dim(x1)[1]
  N=dim(x1)[2]
  mx=matrix(apply(x1,2,mean),1,N)
  Onev=matrix(1,T.,1)
  xrm=x1-Onev%*%mx
  Omega=(1.0/N)*(xrm%*%t(xrm))
 e1=eigen(Omega)
  Ft=e1$vectors[,1:nf]
  SIG=c(1:N)
  for (i in 1:N){
    mm=lm(x1[,i]~Ft)
    sig=sum(mm$residuals^2)/(T.-nf-1)
    SIG[i]=1.0/sqrt(sig)
  }
  DI=diag(SIG)
  y1=x1%*%DI
  ym=matrix(apply(y1,2,mean),1,N)
  y=y1-Onev%*%ym
  Ome=(1.0/N)*(y%*%t(y))
#  
  e2=eigen(Ome)
  Ft=e2$vectors[,1:nf]
  load=matrix(0,N,nf)
  rsq=c(1:N)
  for (i in 1:N){
    mi=lm(x1[,i]~Ft)
    rsq[i]=1-sum(mi$residuals^2)/(var(x1[,i])*(T.-1))
    load[i,]=mi$coefficients[2:(nf+1)]
  }
#print(Ft)
#print(apply(Ft,2,mean))
  cat('Asymptotic PCA:  Extracting', nf, 'factors from',
      T., "observations on", N, 'series\n\n')
  cat('Factor Loading: Summary\n')
#  print('Factor   Minimum  1st-Qu   Median   Mean    3rd-Qu   Maximum')
  lp=matrix(0,nf,7, dimnames=list(paste("F", 1:nf, sep='.'), 
            c("factor", "min", "1st Qu.", "median", "mean", "Q3", "max")))
  for (j in 1:nf){
    lp[j,1]=j
    lp[j, c(2:4, 6:7)] <- quantile(load[, j])
    lp[j, 5] <- mean(load[,j])
  }
  print(lp,digits=4)
  srsq=c(1:6)
  names(srsq) <- dimnames(lp)[[2]][-1] 
  srsq[c(1:3, 5:6)] <- quantile(rsq) 
  srsq[4]=mean(rsq)
  cat('R-squared: Summary\n')
#  print('Min. 1st-Qu  Median  Mean  3rd-Qu  Maximum')
  print(srsq,digits=4)
  apac <- list(eig=e2$values,factors=Ft,loadings=load,rsq=rsq)
  apac
}
