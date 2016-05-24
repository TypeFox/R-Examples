######SCDA ALGORITHM#############
CSDCD.random=function(num.var,num.sample,x.mat,y.vec,d.Star){
  p=num.var
  n=num.sample
  design.Mat=x.mat
  y.sim=y.vec
  lambda=d.Star
  
  ################ DEFINE PARAMETERS ########################
  thru.data=75
  tot.it=thru.data*n
  alpha=matrix(NA,tot.it,n)
  nu.mat=matrix(NA,tot.it,p)
  
  ################ INITIAL VALUES ############################
  
  alpha[1,]=rep(0,n)
  nu.mat[1,]=(alpha[1,]%*%design.Mat*(1/lambda))/n
  
  for(r in 1:thru.data){
    random.seq=sample(1:n, n,replace = FALSE)
    for(k in 1:n){
      t=(r-1)*n+k
      if(r==1){
        t=(r-1)*n+k+1
      }    
      delta.alpha=-(alpha[t-1,random.seq[k]]+t(nu.mat[t-1,])%*%design.Mat[random.seq[k],]-y.sim[random.seq[k]])/(1+(sum(design.Mat[random.seq[k],]**2*(1/lambda))/(2*n)))
      alpha[t,]=alpha[t-1,]+delta.alpha
      nu.mat[t,]=nu.mat[t-1,]+(delta.alpha/n)*design.Mat[random.seq[k],]*(1/lambda)  
    }
  }
  colMeans(nu.mat[(tot.it/2):tot.it,])
}