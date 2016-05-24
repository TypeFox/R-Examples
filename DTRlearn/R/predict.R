predict.linearcl<-function(object,x,...){
  predict=sign(object$bias+x%*%object$beta)
}

predict.rbfcl<-function(object,x,...){
  rbf=rbfdot(sigma=object$sigma)
  n=dim(object$X)[1]
  if (is.matrix(x)) xm=dim(x)[1]
  else if (is.vector(x)) xm=1
  else stop('x must be vector or matrix')
  if (xm==1){ K <- apply(object$X,1,rbf,y=x)
  }else{   K<- matrix(0, xm, n)
    for (i in 1:xm) {K[i,]=apply(object$X,1,rbf,y=x[i,]) }}
  predict=sign(object$bias+K%*%object$alpha1)
}

# predict.rbfcl<-function(object,x){
#   rbf=rbfdot(sigma=object$sigma)
#   n=dim(object$X)[1]
#   if (is.matrix(x)) xm=dim(x)[1]
#   if (is.vector(x)) xm=1
#   K<- matrix(0, xm, n)
#   for (i in 1:n) {
#   K[,i]=apply(x,1,rbf,y=object$X[i,]) }
#   predict=sign(object$bias+K%*%object$alpha1)
# }

predict.qlearn<-function(object,x,...){
  p=dim(x)[2]
  predict=sign(object$co[p+2]+x%*%object$co[(p+3):(2*p+2)])
}
