`hermite` <-
function (x,theta) 
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
  HH<-H(length(theta)-1,x)
  out<-0
  for(k in 1:length(theta)){#k is hermite_k plus 1
  out<-out+theta[k]*HH[[k]]
  }
  out
}

