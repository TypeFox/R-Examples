fwtdr<-function(x,Tij)
{
  
  n = nrow(x)
  d = ncol(x)
  tmp = .Fortran("wtdr",as.matrix(x),as.vector(Tij),as.integer(n),as.integer(d),answer=double(n*d))$answer
  tmp = matrix(tmp,ncol=d)
  return(tmp)
  
  
  
} 
########
### Spatial Rank  Pseudo-Median index
########
spatial_median<-function(rank)
{
  fctindx<-which.min(eucldisc(rank))
  fctindx
}

########
### End Spatial Rank  Pseudo-Median
########

########
### euclidean disitance, give a matirx n*d, return n-vector distance from center
########
eucldisc<-function(x){
  sqx<-x*x
  eucldisc<-apply(sqx,1,sum)
  eucldisc<-sqrt(eucldisc)
  eucldisc
}

########
### End euclidean disitance
########

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

mix_kotz=function(x,mul,sgml,taul){
  #mul is c*d matrix
  #sgml is list
  #taul is vector
  n=dim(x)[1]
  d=dim(x)[2]
  c=dim(mul)[1]   #number of components
  if (sum(taul)!= 1){print("Kotz-EM taul not sum up 1")}
  if (c!=length(sgml)| c!=length(taul)){
    print("Number of component doesn't match!")
    break
  }
  
  Tij=fTij(x,mul,sgml,taul)
    j=0
    param = foreach(j =1:c,.export = c("fwtdr","spatial_median","eucldisc"),.combine="comb",.init=list(list(), list()))%dopar%{ # update median and covariance
  	
  	sgmj = sgml[[j]]
  	eve = eigen(sgmj)
  	evec = eve$vectors
  	evev = eve$values
  	y = x%*%evec%*%diag(evev^(-1/2))%*%t(evec) 	   # standardize
    wtdr=fwtdr(y,Tij[,j])
    ctindx=spatial_median(wtdr)
    mul[j,]= t(x[ctindx,])
       
   ## do loop for covariance, stopping criteria either iter=100 or norm of difference <0.001
   	 s=sum(Tij[,j])    
     x_mius_mu=t(x)-x[ctindx,]  
   
     for (k in 1:100) {
     y_mius_ymu = t(y)-y[ctindx,]
     dis=eucldisc(t(y_mius_ymu))
     invbottom=1/dis
     invbottom[!is.finite(invbottom)]=0    
     sgmj1 = x_mius_mu %*% diag(Tij[,j]/s*invbottom)%*% t(x_mius_mu)+0.00000001*diag(1,d)
     a =norm(sgmj1-sgmj)
     	 if (a<0.0001)
     	 {
       break}
      sgmj = sgmj1
      eve = eigen(sgmj)
  	evec = eve$vectors
  	evev = eve$values
  	y = x%*%evec%*%diag(evev^(-1/2))%*%t(evec)    	
     }  
    sgml[[j]] =sgmj1  
          
   list(mul[j,],sgml[[j]])
  }
  mul<-matrix(unlist(param[[1]]),ncol = d,byrow=T)
  sgml<- param[[2]]
  taul=apply(Tij,2,mean)
  list(mul,sgml,taul)
}

########
### End M-step
########

########
###  em iteration;  criteria (loops<100) or (taul[1] stable)
########

em_kotz_GM=function(x,mul,iter_max){
  numWorkers <- rep(getOption("cl.cores", 2))
  cl <- makeCluster(numWorkers)
  registerDoParallel(cl)
  #x is the data while k is the number of components
  x<-as.matrix(x)
  k<-nrow(mul)
  taull<-list()
  mull<-list()
  sgmll<-list()
  mull[[1]]<-mul
  taull[[1]]<-rep(1/k,k)
  sgmll[[1]]<-lapply(1:k,function(p,q,r) q[[p]]<-diag(nrow=ncol(r)),q=sgmll,r=x)
  
  
  for (t in 1:iter_max){
    mix=mix_kotz(x,mull[[t]],sgmll[[t]],taull[[t]])
    mull=c(mull,list(mix[[1]]))
    sgmll=c(sgmll,list(mix[[2]]))
    taull=c(taull,list(mix[[3]]))
    if (abs(taull[[t+1]][1]-taull[[t]][1])<0.0001)
      {break}
    if(t==iter_max)
    {
      print("The number of iterations reach maxiter and the method may not converge")
    }
  }
  closeAllConnections()
  ## (t+1)=ending position
  list(mull[[t+1]],sgmll[[t+1]],taull[[t+1]],t+1,mull,sgmll,taull)

}