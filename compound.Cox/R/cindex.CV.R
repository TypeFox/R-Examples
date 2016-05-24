cindex.CV <-
function(t.vec,d.vec,X.mat,alpha,K=5){
  p=ncol(X.mat)
  n=length(t.vec)
  PI.cv=numeric(n)
  
  for(k in 1:K){
    temp=(1+(k-1)*n/K):(k*n/K) ### index for the kth fold ###
    
    beta_a_cv=numeric(p)
    for(j in 1:p){
      res=dependCox.reg(t.vec[-temp],d.vec[-temp],X.mat[-temp,j],alpha=alpha,var=FALSE)
      beta_a_cv[j]=res
    }
    
    PI.cv[temp]=X.mat[temp,]%*%beta_a_cv
  }
  
  survConcordance(  Surv(t.vec,d.vec)~PI.cv  )$concordance
}
