

fitclogitboost<-function(Y,M,groupid,iter=100,rho=0.05){
  
  addfunction.var<- vector("list", iter) 
  addfunction<- vector("list", iter)
  theta<- vector("list", iter)
  likeli=c()
  win=Y
  
  
  bigf<-function(x){
    10e-5
  }
  
  for (k in seq(1,iter)){
    
    if (k<2){
      fx=apply(M, 1, bigf)
    }
    
    if (k>1){
      
      bigfa<-function(x){
        temp=0
        for (m in seq(1,k-1)){
          temp=temp+rho*theta[[m]]*predict(addfunction[[m]],x[addfunction.var[[m]]])$y
        }
        return(temp)
      }
      fx=apply(M, 1, bigfa)
    }
    
    
    grad=persamplegrad(fx,win,length(win),groupid)
    
    diff1=rep(1e10,dim(M)[2])
    spllist <- vector("list", dim(M)[2]) 
    predictlist<-vector("list", dim(M)[2]) 
    for(j in seq(1,dim(M)[2])){
      
      spl<- smooth.spline(M[,j],grad)
      spl.predict<-predict(spl,M[,j])$y
      diff1[j]=sum((grad-spl.predict)^2)
      spllist[[j]]=spl
      predictlist[[j]]=spl.predict
    }
    
    addfunction.var[[k]]=which(diff1==min(diff1))
    addfunction[[k]]=spllist[[which(diff1==min(diff1))]]
    addfunction.predict=predictlist[[which(diff1==min(diff1))]]
    
    theta[[k]]=uniclogit(addfunction.predict,win,fx,length(win),groupid)
    
    
    
    likeli=c(likeli,likelihood(fx,win,length(win),groupid))
  }
  
  
  temp=0
  tempm=matrix(0,ncol=dim(M)[2],nrow=length(groupid))
  for (d in seq(1,dim(M)[2])){
    for (j in seq(1,length(groupid))){
      temp=0
      x=M[j,d]
      for (m in seq(1,length(theta))){
        if (addfunction.var[[m]]==d) {
          temp=temp+rho*theta[[m]]*predict(addfunction[[m]],x,deriv=1)$y
        }}
      tempm[j,d]=temp
    }
  }
  
  ederivxj=apply(tempm^2,2,mean)
  varxj=apply(M,2,var)
  rinf=sqrt(ederivxj*varxj)
  
  mmax=apply(M,2,max)
  mmin=apply(M,2,min)
  
  return(list(func=addfunction,index=addfunction.var,theta=theta,likeli=likeli,rinf=rinf,rho=rho,xmax=mmax,xmin=mmin))
}
