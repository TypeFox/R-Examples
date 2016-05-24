dqreg <-
function(object,y,X){
  n<-length(y)
  X<-matrix(X,n,length(X),byrow=T)
  beta<- apply(object$beta,2,median)
  alpha<-apply(object$alpha,2,median)
  alpha<-t(matrix(alpha,nrow=object$L))
  shape<-median(object$shape)
  d<-log_pdf(y,X,beta,alpha,object$kappa,object$base,shape)
  return(exp(d))}
