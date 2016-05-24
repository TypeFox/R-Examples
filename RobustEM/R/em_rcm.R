#########################################################################################################
#The spatial-EM algorithm. This algorithm accepts as input the data x and returns as output the
#mean and covariance of each component of a distribution. The mean is calculated using the spatial
#median and it is represented as mul, while the covariance is calculated using the Modified Rank 
#Covariance Matrix (MRCM)
#########################################################################################################

fwtdr<-function(x,Tij)
{
  
  n = nrow(x)
  d = ncol(x)
  tmp = .Fortran("wtdr",as.matrix(x),as.vector(Tij),as.integer(n),as.integer(d),answer=double(n*d))$answer
  tmp = matrix(tmp,ncol=d)
  return(tmp)
  
} 

spatial_median<-function(rank)
{
   fctindx<-which.min(eucldisc(rank))
   fctindx
}

# euclidean distance, give a matrix n*d, return n-vector distance from center

eucldisc<-function(x){
  sqx<-x*x
  eucldisc<-apply(sqx,1,sum)
  eucldisc<-sqrt(eucldisc)
  eucldisc<-as.vector(eucldisc)
  eucldisc
}

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
### M-Step. The mix_rcm uses the value gotten in the E-step to compute the spatial median and the rank covariance
#########
#################################################################################
#function to ensure that for each can print more than one variable
################################################################################
comb <- function(x, ...) {
  lapply(seq_along(x),function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

mix_rcm<-function(x,mul,sgml,taul){
  #mul is c*d matrix, sgml is list, taul is vector
  n=nrow(x)
  d=ncol(x)
  c=nrow(mul)   #number of components
 
  if (c!=length(sgml)| c!=length(taul)){
    stop("Number of component doesn't match!")
    
  }
  Tij<-fTij(x,mul,sgml,taul)
  
    j=0
    param<-foreach(j = 1:c,.export = c("fwtdr","spatial_median","eucldisc"),.combine="comb",.init=list(list(), list()))%dopar%{
    wtdr<-fwtdr(x,Tij[,j])
    ctindx <- spatial_median(wtdr)
    mul[j,]<-x[ctindx,]   #mean
    
    ## Use weighted MRCM find Sigma (eigenvalue=mean abs deviation)
    s<-sum(Tij[,j])
    #R<-wtdr[[1]]
    wtdrcm<-t(diag(Tij[,j]/s)%*%wtdr)%*%wtdr#covariance at t+1
    evec<-eigen(wtdrcm)$vectors
    lambda<-rep(0,d)
    
    ##Use Kai's weighted lambda
    jthcomp_indx<-ceiling(n*(1-s/n))
    
    for (k in 1:d){
      val<-t(mul[j,])%*%evec[,k]
      
      wtdprj<-Tij[,j]*(x%*%evec[,k]-rep(val,n))
      abswtdprj<-abs(wtdprj)
      wtdprj_jthcomp<-wtdprj[which(abswtdprj>=sort(abswtdprj)[jthcomp_indx])]
      lambda[k]<-median(abs(wtdprj_jthcomp))*1.4826 #1.4826 consistent coef.
      
    }
    
    ##End Use Kai's weightd lambda
    
    mrcm<-evec%*%diag(lambda^2)%*%t(evec)
    
    
    sgml[[j]]<-mrcm+0.00000001*diag(1,d)
    
    list(mul[j,],sgml[[j]])
    
  }
  mul<-matrix(unlist(param[[1]]),ncol = d,byrow=T)
  sgml<- param[[2]]
  taul<-apply(Tij,2,mean)
  list(mul,sgml,taul)
}
########
### End M-step
########


em_rcm_GM<-function(x,mul,iter_max){
  
  numWorkers <- rep(getOption("cl.cores", 2))
  cl <- makeCluster(numWorkers)
  registerDoParallel(cl)
  
  x<-as.matrix(x)
  k<-nrow(mul)
  taull<-list()
  mull<-list()
  sgmll<-list()
  mull[[1]]<-mul
  taul = rep(1/k,k)
  taull[[1]]<-taul
  sgmll[[1]]<-lapply(1:k,function(p,q,r) q[[p]]<-diag(nrow=ncol(r)),q=sgmll,r=x)
  
  
  
  
  for (t in 1:iter_max){
    mix<-mix_rcm(x,mull[[t]],sgmll[[t]],taull[[t]])
    mull<-c(mull,list(mix[[1]]))
    sgmll<-c(sgmll,list(mix[[2]]))
    taull<-c(taull,list(mix[[3]]))
    if (abs(taull[[t+1]][1]-taull[[t]][1])<0.0001)
    {break}
    
    if(t==iter_max)
    {
      print("The number of iterations reach maxiter and the method may not converge")
    }
    
  }
  closeAllConnections()
  ## (t+1)=ending position
  return(list(mull[[t+1]],sgmll[[t+1]],taull[[t+1]],t+1,mull,sgmll,taull))
}

