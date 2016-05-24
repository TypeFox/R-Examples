Plearning<-function(X,AA,RR,n,K,pi,pentype='lasso',kernel='linear',sigma=c(0.03,0.05,0.07),clinear=2.^(-2:2),m=4,e=0.00001){
select=matrix(1,n,1)
QL=matrix(0,n,K)
M=matrix(1,n,K)
C=matrix(1,n,K)
models=list()
prob=matrix(1,n,K)
QLproj=matrix(0,n,K+1)
Qspecify=matrix(0,n,K)
QR_future=0
Rsum=0
if (is.matrix(X)){
for (k in K:1){
A=AA[[k]]
output_Q=Qlearning_Single(X,A,RR[[k]]+QR_future,pentype=pentype,m=m)
QR_future=output_Q$Q
#subsititute the outcome by expected outcome of best treatment
QL[,k]=output_Q$Q
if(k<K) R_p=Rsum*select/prob[,K]+apply(QLproj[,(k+1):K]%*%Qspecify[,(k+1):K],2,sum)
if(k==K) R_p=Rsum*select/prob[,K]
R=(RR[[k]]+R_p)
if (kernel=='linear'){
models[[k]]=Olearning_Single(X,A,R,pi[[k]],pentype=pentype,clinear=clinear,e=e,m=m)
}else if (kernel=='rbf'){
  models[[k]]=Olearning_Single(X,A,R,pi[[k]],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
}else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))

right=(sign(models[[k]]$fit)==A)
#update fo next stage
M[,k:K]=M[,k:K]*(right%*%rep(1,K-k+1))
if (k>1) C[,k:K]=M[,k-1:K-1]-M[,k:K]
if (k==1){
  C[,2:K]=M[,1:(K-1)]-M[,2:K]
  C[,1]=rep(1,n)-M[,1]
}

select=select*right
prob[,k:K]=prob[,k:K]*(pi[[k]]*rep(1,K-k+1))
Rsum=rep(1,n)
for (j in k:K){
if (j>1) {QLproj[,j]=(C[,j]-(1-pi[[j]])*M[,j-1])/prob[,j]
} else QLproj[,1]=(C[,j]-(1-pi[[j]]))/prob[,j]

Qspecify[,j]=QL[,j]+Rsum
Rsum=Rsum+RR[[j]]
}}}

if (is.list(X)){
  for (k in K:1){
    A=AA[[k]]
    output_Q=Qlearning_Single(X[[k]],A,RR[[k]]+QR_future,pentype=pentype)
    QR_future=output_Q$Q
    #subsititute the outcome by expected outcome of best treatment
    QL[,k]=output_Q$Q
    if(k<K) R_p=Rsum*select/prob[,K]+apply(QLproj[,(k+1):K]%*%Qspecify[,(k+1):K],2,sum)
    if(k==K) R_p=Rsum*select/prob[,K]
    R=(RR[[k]]+R_p)
    if (kernel=='linear'){
      models[[k]]=Olearning_Single(X[[k]],A,R,pi[[k]],pentype=pentype)
    }else if (kernel=='rbf'){
      models[[k]]=Olearning_Single(X[[k]],A,R,pi[[k]],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
    }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))
    
    right=(sign(models[[k]]$fit)==A)
    #update fo next stage
    M[,k:K]=M[,k:K]*(right%*%rep(1,K-k+1))
    if (k>1) C[,k:K]=M[,k-1:K-1]-M[,k:K]
    if (k==1){
      C[,2:K]=M[,1:(K-1)]-M[,2:K]
      C[,1]=rep(1,n)-M[,1]
    }
    
    select=select*right
    prob[,k:K]=prob[,k:K]*(pi[[k]]*rep(1,K-k+1))
    Rsum=rep(1,n)
    for (j in k:K){
      if (j>1) {QLproj[,j]=(C[,j]-(1-pi[[j]])*M[,j-1])/prob[,j]
      } else QLproj[,1]=(C[,j]-(1-pi[[j]]))/prob[,j]
      
      Qspecify[,j]=QL[,j]+Rsum
      Rsum=Rsum+RR[[j]]
    }}}
models
}