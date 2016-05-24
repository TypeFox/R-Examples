delta_k <-
function(Y,X,D,lambda1,lambda2){
   gamma<-0.5
   2*sum(as.vector(X-Y)%*%as.vector(D))+2*gamma*sum(D^2)+
     P(X+D,lambda1,lambda2)-P(X,lambda1,lambda2)
}