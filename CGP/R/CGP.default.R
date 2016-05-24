CGP.default <-
function(X,yobs,...){
  
  X<-as.matrix(X)
  yobs<-as.numeric(yobs)
  
  est<-CGPEst(X,yobs)
  est$call<-match.call()
  class(est)<-"CGP"
  return(est)
    
}
