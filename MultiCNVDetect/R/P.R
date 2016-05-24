P <-
function(X,lambda1,lambda2){
  p<-dim(X)[2]
  return(lambda1*sum(apply(X,2,norm_2))+
    lambda2*sum(apply(X[,1:(p-1)]-X[,2:p],2,norm_2)))
}