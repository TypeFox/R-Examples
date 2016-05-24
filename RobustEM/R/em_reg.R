#########################################################################################################
#The regular-EM algorithm. This algorithm accepts as input the data x and returns as output the
#mean and covariance of each component of a distribution.
#########################################################################################################

#########
### E-step, return  P(Z_i=j|X,theta) as n vector  for Gussian Dist.
#########
fTij<-function(x,mul,sgml,taul)
{
  c = nrow(mul)
  fTij<-matrix(0,nrow=nrow(x),ncol=c)
  fTij<-sapply(1:c,function(z,x,m,s) dmvnorm(x,m[z, ],s[[z]]),x=x,m=mul,s=sgml)
  fTij[fTij==0]=1
  fTij = fTij%*%diag(taul)
  denom<-apply(fTij,1,sum)
  fTij<-apply(fTij,2,function(x) x/denom)
  fTij
  
}

#########
### End E-step
#########

#########
### M-Step
#########

mix_reg=function(x,mul,sgml,taul){
  #mul is c*d matrix, sgml is list, taul is vector
  n=nrow(x)
  d=ncol(x)
  c=nrow(mul)   #number of components
#   if (sum(taul)!= 1)
#   {print("reg EM taul not sum up 1")}
  if (c!=length(sgml)| c!=length(taul)){
    stop("Number of component doesn't match!")
    
  }
  Tij=fTij(x,mul,sgml,taul)
  

  for (j in 1:c){    # update mean and covariance
    s=sum(Tij[,j])
    mul[j,]=t(x)%*%Tij[,j]/s
    mulj=as.vector(mul[j,])

    x_minus_mu=t(x)-mulj
    sgml[[j]]=x_minus_mu%*%diag(Tij[,j]/s)%*%t(x_minus_mu)+0.00000001*diag(1,d)
  }
  taul=apply(Tij,2,mean)
  list(mul,sgml,taul)
}

########
### End M-step
########

########
###  em iteration;  criteria (loops<100) & (taul[1] stable)
########

em_reg_GM=function(x,mul,iter_max){
  x<-as.matrix(x)
  k<-nrow(mul)
  taull<-list()
  mull<-list()
  sgmll<-list()
  mull[[1]]<-mul
  taull[[1]]<-rep(1/k,k)
  sgmll[[1]]<-lapply(1:k,function(p,q,r) q[[p]]<-diag(nrow=ncol(r)),q=sgmll,r=x)
  
  for (t in 1:iter_max){
    mix=mix_reg(x,mull[[t]],sgmll[[t]],taull[[t]])
    mull=c(mull,list(mix[[1]]))
    sgmll=c(sgmll,list(mix[[2]]))
    taull=c(taull,list(mix[[3]]))
 
    if (abs(taull[[t+1]][1]-taull[[t]][1])<0.0001)
      {break}
  }
  ## (t+1)=ending position
  if(t==iter_max)
  {
    print("The number of iterations reach maxiter and the method may not converge")
  }
  list(mull[[t+1]],sgmll[[t+1]],taull[[t+1]],t+1,mull,sgmll,taull)

}
