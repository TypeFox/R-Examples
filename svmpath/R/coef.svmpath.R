coef.svmpath<-function(object,lambda,...){
  if(missing(lambda)){
    alpha<-object$alpha
    lambda<-object$lambda
    alpha0<-object$alpha0
    }
  else{
    alphs<-predict(object,lambda=lambda,type="alpha")
    alpha<-alphs$alpha
    alpha0<-alphs$alpha0
    }
  alpha<-alpha*object$y
  beta<-scale(t(object$x)%*%alpha,FALSE,lambda)
  beta0<-alpha0/lambda
  list(beta=beta,beta0=beta0,lambda=lambda)
}
  
