choosePriors <-
function(kmax,priors){
	index=which(priors[,1]==kmax)

	if(length(index)==0){	
		print("Sorry, distributions are avaiable for kmax (maximum number of Changepoints/Edges) between 2 and 40 only.")
	}
			
	if(length(index)>1){
		par(mfrow=(c(ceiling(length(index)/2),2)))
		for(i in index[2:1]){
			plot(0:kmax,priors[i,4:(kmax+4)],type="h",lwd=5,col=2,ylim=c(0,1),main=paste("alpha=",priors[i,2],", beta= ", priors[i,3]),ylab="Prior probability",xlab="Number of changepoints or TF")
		}
		for(i in index[3:length(index)]){
			plot(0:kmax,priors[i,4:(kmax+4)],type="h",lwd=5,col=4,ylim=c(0,1),main=paste("alpha=",priors[i,2],", beta= ", priors[i,3]),ylab="Prior probability",xlab="Number of changepoints or TF")
		}
	}

	if(length(index)==1){	
		par(mfrow=(c(1,1)),cex=1)
		plot(0:kmax,priors[index,4:(kmax+4)],type="h",lwd=5,col=2,ylim=c(0,1),main=paste("alpha=",priors[index,2],", beta= ", priors[index,3]),ylab="Prior probability",xlab="Number of changepoints or TF")
	}
	
}
