Olearning_Single<-function(H,A,R2,pi=rep(1,n),pentype='lasso',kernel='linear',sigma=c(0.03,0.05,0.07),clinear=2.^(-2:2),m=4,e=1e-5)
{
npar=length(clinear)
n=length(A)
p=dim(H)[2]

if (max(R2)!=min(R2)){
  if (pentype=='lasso'){
    cvfit=cv.glmnet(H,R2,nfolds=m)
    co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
  }else if (pentype=='LSE'){
    co=coef(lm(R2~H))
  }else stop(gettextf("'pentype' is the penalization type for the regression step of Olearning, the default is 'lasso',
it can also be 'LSE' without penalization"))
    r=R2-cbind(rep(1,n),H)%*%co
  } else r=R2
rand=sample(m,n,replace=TRUE)

r=r/pi
if (kernel=='linear'){
V=matrix(0,m,npar)
for (i in 1:m){
this=(rand!=i)
X=H[this,]
Y=A[this]
R=r[this]
Xt=H[!this,]
Yt=A[!this]
Rt=r[!this]
for (j in 1:npar){
model=wsvm(X,Y,R,C=clinear[j],e=e)
YP=predict(model,Xt)
V[i,j]=sum(Rt*(YP==Yt))/sum(YP==Yt)
}}
mimi=colMeans(V)
best=which.max(mimi)
cbest=clinear[best]
model=wsvm(H,A,r,C=cbest,e=e)}

if (kernel=='rbf'){
  nsig=length(sigma)
  V=array(0,c(npar,nsig,m))
  for (i in 1:m){
  this=(rand!=i)
  X=H[this,]
  Y=A[this]
  R=r[this]
  Xt=H[!this,]
  Yt=A[!this]
  Rt=r[!this]
for (j in 1:npar){
  for (s in 1:nsig){
    model=wsvm(X,Y,R,'rbf',sigma=sigma[s],C=clinear[j],e=e)
    YP=predict(model,Xt)
    V[j,s,i]=sum(Rt*(YP==Yt))/sum(YP==Yt)
  }}}
mimi=apply(V,c(1,2),mean)
best=which(mimi==max(mimi),arr.ind=TRUE)
bestC=clinear[best[1]]
bestSig=sigma[best[2]]
print(bestC)
print(bestSig)
model=wsvm(H,A,r,'rbf',bestSig,C=bestC,e=e)}
model
}

Olearning<-function(X,AA,RR,n,K,pi,pentype='lasso',kernel='linear',sigma=c(0.03,0.05,0.07),clinear=2.^(-2:2),m=4,e=0.00001){
select=rep(TRUE,n)
R_future=0
prob=rep(1,n)
models=list()
if (is.matrix(X)){
for (j in K:1){
  R=(RR[[j]]+R_future)/(prob*pi[[j]])  
  if (kernel=='linear'){
    models[[j]]=Olearning_Single(X[select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,clinear=clinear,e=e,m=m)
  }else if (kernel=='rbf'){
    models[[j]]=Olearning_Single(X[select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
  }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))
  
  
  
  select[which(select==1)]=(sign(models[[j]]$fit)==AA[[j]][select])
  R_future=R_future+RR[[j]]
  prob=prob*pi[[j]]
}
}

if (is.list(X)){
  for (j in K:1){
    R=(RR[[j]]+R_future)/(prob*pi[[j]])
    if (kernel=='linear'){
      models[[j]]=Olearning_Single(X[[j]][select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,clinear=clinear,e=e,m=m)
    }else if (kernel=='rbf'){
      models[[j]]=Olearning_Single(X[[j]][select,],AA[[j]][select],R[select],pi[[j]][select],pentype=pentype,kernel=kernel,sigma=sigma,clinear=clinear,e=e,m=m)
    }else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))      
    select[which(select==1)]=(sign(models[[j]]$fit)==AA[[j]][select])
    R_future=R_future+RR[[j]]
    prob=prob*pi[[j]]
  }
}

class(models)=paste('Olearning',pentype,sep='')
# sumR=RR[[1]]
# if (K>=2){
# for (j in 2:K) sumR=sumR+RR[[j]]
# }
# value=sum(sumR*select)/sum(select)
models
}