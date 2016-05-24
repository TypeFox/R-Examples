Simu.PDMP <-
function(x0,T,verbose=TRUE){
	S<-c(); x<-x0
	for (i in 1:T){
		a<-Simu.Cond.HR(1,.JumpRate,x[i],verbose=FALSE)
		s<-min(a,.tstar(x[i]))
		S<-c(S , s)
		x<-c(x,.Transition(x[i],S[i]))
	}
	Retour<-matrix(0,nrow=T+1,ncol=2); Retour[,1]<-x; Retour[1:T,2]<-S
	if (verbose){
		plot(c(0,cumsum(S)),x,"s",main="Simulations",xlab="Time",ylab="State space")
	}
	Retour
}
