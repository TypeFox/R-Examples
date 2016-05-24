Simu.Cond.HR <-
function(N , lambda , x , verbose=TRUE){
	rep<-c()
	for (k in 1:N){
		a<-1/20; c<-0; A<-0
		while (c==0){
			M<-max( lambda( x , seq(from=A,to=A+a,by=a/10) ) )
			n<-rpois(1 , lambda=M*a)
			if (n>0){
				d<-0
				while (d==0){
					u<-runif(1,min=A,max=A+a)
					v<-runif(1,min=0,max=M)
					if (v<lambda( x , u)){
						rep<-c(rep,u); d<-1
					}
				}
			c<-1
			} else{A<-A+a}
		}
	}
	if (verbose){
		if ( length(x) > 1){
			hist(rep,probability=1,main=paste("Histogram of simulations"),xlab="Simulations",ylab=paste("Conditional density given X=(",paste(x,collapse=","),")",sep=""));
			} else {hist(rep,probability=1,main=paste("Histogram of simulations"),xlab="Simulations",ylab=paste("Conditional density given X=",x,sep=""));}
		print(paste("mean =",mean(rep))); print(paste("variance =",var(rep)));
		}
	rep
}
