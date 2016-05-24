`hermite.coef` <-
function (z,kk) 
{
  H<-function(kk,y){
  HH<-list(kk+1)
  HH[[0+1]]<-rep(1,length(y))
  HH[[1+1]]<--y
  if(kk>1){
  for(k in 2:kk){
  HH[[k+1]]<--1/sqrt(k)*y*HH[[k+1-1]]-sqrt((k-1)/k)*HH[[k+1-2]]
  }
  }
  HH
  }
  len<-length(z)
  z<-jitter(z)
  z<-sort(z)
  q<-table(z)
  q<-table(z)/len
  theta<-numeric(kk+1)
  theta[0+1]<-sum(q*z)
  gy<-seq(0,length=len, by=1/len)
  gy<-qnorm(gy)
  y<-dnorm(gy)
  y[c(1,len)]<-0
  HH<-H(kk,gy)
  for(k in 1:kk){
  theta[k+1]<-sum(
  (z[1:(len-1)]-z[2:len])/sqrt(k)*HH[[k+1-1]][-1]*y[-1]
  )
  }
  theta
}

