Zscore.internal <-
function(Theta, Tau, X){
	
	if(is.vector(X)) X<- t(as.matrix(X))
	
	N<- nrow(X)
	M<- ncol(X)
	G<- length(Tau)
	
	dum<-array(apply(Theta,1,dbinom, size=1, x=t(X)), dim=c(M,N,G))
	
	Z1<-t(Tau*t(apply(dum, c(2,3), prod)))
	Z<-Z1/apply(Z1,1,sum)
	
	Z
	
	}
