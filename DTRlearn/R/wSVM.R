wsvm<-function(X,A,wR,kernel='linear',sigma=0.05,C=1,e=0.00001){
  wAR=A*wR
  if (kernel=='linear'){
    K=X%*%t(X)
  }
  else if (kernel=='rbf'){
    rbf=rbfdot(sigma=sigma)
    K=kernelMatrix(rbf,X)
  } else stop(gettextf("Kernel function should be 'linear' or 'rbf'"))
  
  H=K*(wAR%*%t(wAR))
  n=length(A)
  solution=ipop(-abs(wR),H,t(A*wR),0,numeric(n),rep(C,n),0,maxiter=100)
  alpha=primal(solution)
alpha1=alpha*wR*A
if (kernel=='linear'){
w=t(X)%*%alpha1 #parameter for linear
fitted=X%*%w
rm=sign(wR)*A-fitted
} else if (kernel=='rbf'){
  #there is no coefficient estimates for gaussian kernel
  #but there is fitted value, first we compute the fitted value without adjusting for bias
  fitted=K%*%alpha1 
  rm=sign(wR)*A-fitted
} 
Imid =(alpha < C-e) & (alpha > e)
rmid=rm[Imid==1];
if (sum(Imid)>0){
  bias=mean(rmid)
} else {
    Iup=((alpha<e)&(A==-sign(wR)))|((alpha>C-e)&(A==sign(wR)))
    Ilow=((alpha<e)&(A==sign(wR)))|((alpha>C-e)&(A==-sign(wR)))
    rup=rm[Iup]
    rlow=rm[Ilow]
    bias=(min(rup)+max(rlow))/2}
fit=bias+fitted
if (kernel=='linear') {
  model=list(alpha1=alpha1,bias=bias,fit=fit,beta=w)
  class(model)<-'linearcl'
} else if (kernel=='rbf') {
  model=list(alpha1=alpha1,bias=bias,fit=fit,sigma=sigma,X=X)
  class(model)<-'rbfcl'}
return (model)
}

