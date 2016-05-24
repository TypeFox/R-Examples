Fx <-
function(Y,X,lambda1,lambda2){
  return(sum(as.vector((Y-X))^2)+P(X,lambda1,lambda2))
}