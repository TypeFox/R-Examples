estim_herit <- function(Y,Z)
{
  n=length(Y)
  N=ncol(Z)
  a=n/N
  nb_iter=20
  eta_init=c(0.1,0.5,0.9)
  eta_mat=matrix(0,3,nb_iter+1)
  for (i in 1:3)
  {
    eta_mat[i,1]=eta_init[i]
  }
  #Z=scale(W,center=TRUE,scale=TRUE)
  #M=Z%*%t(Z)/N
  M=prod_cpp(Z)/N
  eig_M=eigen(M)
  O=eig_M$vectors
  lambda=eig_M$values
  Y_tilde=t(O)%*%Y
  for (i in 1:3)
  {
    for (nb in 1:nb_iter)
    {
      A=sum(Y_tilde^2*(lambda-1)/(eta_mat[i,nb]*(lambda-1)+1)^2)/sum(Y_tilde^2/(eta_mat[i,nb]*(lambda-1)+1))-1/n*sum((lambda-1)/(eta_mat[i,nb]*(lambda-1)+1))
      B=((-2*sum(Y_tilde^2*(lambda-1)^2/(eta_mat[i,nb]*(lambda-1)+1)^3)*sum(Y_tilde^2/(eta_mat[i,nb]*(lambda-1)+1))+(sum(Y_tilde^2*(lambda-1)/(eta_mat[i,nb]*(lambda-1)+1)^2))^2)/(sum(Y_tilde^2/(eta_mat[i,nb]*(lambda-1)+1)))^2 +1/n*sum((lambda-1)^2/(eta_mat[i,nb]*(lambda-1)+1)^2))
      eta_mat[i,(nb+1)]=eta_mat[i,nb]-A/B
    }
  }
  
  ind_init=which.max(apply(cbind(abs(eta_mat[,(nb+1)]-0.01),abs(eta_mat[,(nb+1)]-0.99)),1,min))
  
  eta_chap=eta_mat[ind_init,(nb+1)]
  
  sigma2_chap=max(1/n*sum(Y_tilde^2/(eta_chap*(lambda-1)+1)),10^-6)
  list(heritability=min(max(eta_chap,0),1),sig2=sigma2_chap)
}
