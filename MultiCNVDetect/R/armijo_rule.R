armijo_rule <-
function(Y,X,D,lambda1,lambda2,alpha){
  gamma<-0
  beta<-0.5
  sigma<-0.1
  temp<-min(alpha/beta,1)
  alpha<-temp
  
  j<-1
  a<-Fx(Y,X,lambda1,lambda2)
  b<-sigma*delta_k(Y,X,D,lambda1,lambda2)
  while( Fx(Y,X+alpha*D,lambda1,lambda2)> a+alpha*b ){
    alpha<-temp*(beta^j)
    j<-j+1
    }
  return(alpha)
}