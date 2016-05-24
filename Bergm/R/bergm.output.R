bergm.output <- function(x,
                         ...){
 
if(x$nchains > 1){ # ADS == TRUE
	
	cat("\n", "MCMC results for Model: y ~", paste(x$formula[3]),"\n")
	
	cat("\n","Posterior mean:","\n")
	postmean <- apply(x$Theta,c(3,2),mean)
	rownames(postmean) <- paste("Chain",seq(1,x$nchains)," ")
	colnames(postmean) <- paste("theta",seq(1,x$dim)," (",
	                            x$specs[seq(1,x$dim)],")",sep="")
	postmean <- as.table(postmean)
	print(postmean)

	cat("\n",'Posterior sd:',"\n")
	postmean <- apply(x$Theta,c(3,2),sd)
	rownames(postmean) <- paste("Chain",seq(1,x$nchains)," ")
	colnames(postmean) <- paste("theta",seq(1,x$dim)," (",
	                            x$specs[seq(1,x$dim)],")",sep="")
	postmean <- as.table(postmean)
	print(postmean)

	cat("\n")
	rates <- matrix(x$acc.rate,x$nchains,1)
	rownames(rates) <- paste("Chain",seq(1,x$nchains)," ")
	colnames(rates) <- paste("Acceptance rate:")
	rates <- as.table(rates)
	print(rates)
	cat("\n")
	
	FF <- apply(x$Theta,2,cbind)
	
	cat("\n",'Overall posterior density estimate:',"\n")
	overall <- rbind(apply(FF,2,mean),apply(FF,2,sd))
	rownames(overall) <- c("Post. mean","Post. sd")
	colnames(overall) <- paste("theta",seq(1,x$dim)," (",
	                           x$specs[seq(1,x$dim)],")",sep="")
	all <- as.table(overall)
	print(overall)
		
}else{ # ADS == FALSE
	
	FF <- x$Theta
	
	cat("\n", "MCMC results for Model: y ~", paste(x$formula[3]),"\n")
	
	cat("\n",'Posterior density estimate:',"\n")
	overall <- rbind(apply(FF,2,mean),apply(FF,2,sd))
	rownames(overall) <- c("Post. mean:","Post. sd:")
	colnames(overall) <- paste("theta",seq(1,x$dim)," (",
	                           x$specs[seq(1,x$dim)],")",sep="")
	all <- as.table(overall)
	print(overall)
	
	rates <- matrix(x$acc.rate,1,x$dim)	

}		

	dev.new()
	
	seqq <- 4 # maximum number of parameter diagnostics plots to display in each graphics window
    	
	K <- mcmc(data=FF)
	par(mfrow=c(min(4,x$dim),3),oma=c(0,0,3,0),mar=c(4,3,1.5,1))

	for(i in 1:x$dim){
		if(i%in%c(5,9,13)){
			dev.new()
			par(mfrow=c(min(4,x$dim-(i-1)),3),
			            oma=c(0,0,3,0),
			            mar=c(4,3,1.5,1))
		}
		plot(density(FF[,i]),
		     main="",
		     axes=FALSE, 
		     xlab=bquote(paste(theta[.(i)]," (",.(x$specs[i]),")")),
		     ylab="",lwd=2)
		axis(1); axis(2)
		plot(FF[,i],type="l",xlab="Iterations",ylab="")
		autocorr.plot(K[,i],auto.layout=FALSE,...)
		
		if(x$dim > 4) 
			seqq <- seq(4, x$dim, 4)
		if(i%in%union(x$dim,seqq)) 
			title(paste("MCMC output for Model: y ~",x$formula[3]),outer=TRUE)
	
}
cat("\n","Overall acceptance rate:", format(mean(rates), digits = 2),"\n","\n","\n")
}

