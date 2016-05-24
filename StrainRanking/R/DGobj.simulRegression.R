

## constructor of DGobj based on a mechanistic model
DGobj.simul.regression=function(sqrtn, size1, size2, theta, alpha.function, sigma, plots=FALSE){
	
	## Definition of the spatial grid
	N=sqrtn^2
	x0 <- seq(1,sqrtn,length=sqrtn)
	x <- expand.grid(x0,x0)

	## Simulation of parameters of the Dirichlet distribution
	nbStrains=length(theta)
	alpha <- alpha.function(x)
	alpha0 <- rowSums(alpha)
	
	if(ncol(alpha)!=length(theta)){
		stop("[simul.regression] Error: The length of theta should be equal to the number of columns produced by alpha.function()")
	} else {}

	## Simulation of strain proportions from the Dirichlet distribution
	pS <- NULL
	for (i in 1:nrow(alpha)) {
	  	pS <- rbind(pS,as.numeric(.rdirichlet(1,alpha=alpha[i,])))
	}

	## Simulation of the random noise
	eta <- rnorm(N,0,sigma)

	## Computation of the growth variables
	Z <- eta
	for (j in 1:ncol(pS)) {
	  	Z <- Z + (theta[j] * pS[,j])
	}

	## Sampling
	echantillon <- sample(1:N,size1)
	tirage <- NULL
	for (i in echantillon) {
	  tirage <- cbind(tirage,rmultinom(1,size2,pS[i,]))
	}

	if(plots){
		## Production of proportion plots
		par(mfrow=c(1,nbStrains),mar=c(2,2,1,0.2))
		for(j in 1:nbStrains){
			plot(x,cex=pS[,j]*3,pch=19,asp=1,axes=F,xlab="",ylab="",main=paste("Prop. strain ",j))
		}
	}
	
	## Construction of the genetic and demographic data sets
	geneticData=cbind(x[echantillon,],t(tirage))
	demographicData=cbind(x,Z)
	colnames(geneticData)=c("x1","x2",paste("strain",1:nbStrains,sep=""))
	colnames(demographicData)=c("x1","x2","growth")
	
    return(new("DGobj",demographic=as.matrix(demographicData),genetic=as.matrix(geneticData)))

}






