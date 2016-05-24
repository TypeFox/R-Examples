Qlearning_Single<-function(H,A,R,pentype='lasso',m=4){
n=length(A)
X=cbind(H,A,diag(A)%*%H)
if (pentype=='lasso'){
  cvfit=cv.glmnet(X,R,nfolds=m)
  co=as.matrix(predict(cvfit,s="lambda.min",type="coeff"))
}else if (pentype=='LSE'){
  co=coef(lm(R~X))
}else stop(gettextf("'pentype' is the penalization type for the regression step of Olearning, the default is 'lasso',
it can also be 'LSE' without penalization"))

XX1=cbind(rep(1,n),H,rep(1,n),diag(n)%*%H)
XX2=cbind(rep(1,n),H,rep(-1,n),-diag(n)%*%H)
Q1=XX1%*%co
Q2=XX2%*%co
Q=apply(cbind(Q1,Q2),1,max)
Qsingle=list(co=co,Q=Q)
class(Qsingle)='qlearn'
Qsingle
}

Qlearning<-function (X,AA,RR,K,pentype='lasso',m=4) {
R_future=0
coef=list()
models=list()
if (is.matrix(X)){
for (j in K:1){
R=RR[[j]]+R_future
if (min(R)!=max(R)){
models[[j]]=Qlearning_Single(X,AA[[j]],R,pentype=pentype)
R_future=models[[j]]$Q}
else {
  models[[j]]=list(co=rep(0,2+2*dim(X)[2]),Q=R)
  R_future=R
}}
}

if (is.list(X)){
  for (j in K:1){
    R=RR[[j]]+R_future
    if (min(R)!=max(R)){
      models[[j]]=Qlearning_Single(X[[j]],AA[[j]],R,pentype=pentype,m=4)
      R_future=models[[j]]$Q}
    else {
      models[[j]]=list(co=rep(0,2+2*dim(X[[j]])[2]),Q=R)
      R_future=R }}
}
models}