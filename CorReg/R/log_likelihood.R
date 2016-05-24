# ' logvraisemblance d'une sous-r?gression
# '
# ' @param theta vector of the parameters c(sigma,beta)
log_likelihood<-function(theta=theta,Y=Y,X=X){
  X=as.matrix(X)#cbind(rep(1,times=nrow(X)),X)
  return(-(length(Y)/2)*log(2*pi*theta[1]^2)-(sum((Y-as.matrix(X)%*%theta[-1])^2))/(2*theta[1]^2) ) 
}