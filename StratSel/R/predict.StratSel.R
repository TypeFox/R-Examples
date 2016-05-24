predict.StratSel <-
function(object, prob=FALSE, profile, ...){
 	
 	x <- object
 	if(missing(profile)){
 		dim.x11 <- x$DIM[c(1:2)]
 		dim.x14 <- x$DIM[c(3:4)]
 		dim.x24 <- x$DIM[c(5:6)]
 	
		X11 <- model.matrix(x$f, data = x$model, rhs = 1)
		X14 <- model.matrix(x$f, data = x$model, rhs = 2)
		X24 <- model.matrix(x$f, data = x$model, rhs = 3)
	
	 	beta <- x$coefficient
	 	# predicting outcomes
		a <- dim.x11[2] + dim.x14[2] + 1
		b <- length(beta)-1
		if (x$corr==FALSE) b <- length(beta)

		y24.latent <- X24%*%as.matrix(beta)[a:b]
		p24 <- pnorm(y24.latent)
		a <- dim.x11[2]+1
		b <- dim.x11[2]+dim.x14[2]
		y14.latent <- p24*X14%*%as.matrix(beta)[a:b] - X11%*%as.matrix(beta)[1:dim.x11[2]]
	
		rho <- beta[length(beta)]

###
		part1 <- cbind(-y24.latent,y14.latent); part2 <- cbind(y24.latent,y14.latent)		# combining for mvnorm
		covar1 <- matrix(c(1,-rho,-rho,1),2,2)				# negative rho
		covar2 <- matrix(c(1,rho,rho,1),2,2)	
		covar1.nc <- matrix(c(1,0,0,1),2,2)				# negative rho
		covar2.nc <- matrix(c(1,0,0,1),2,2)	
	
		prob3 <- 1-pnorm(y14.latent)	#  outcome 1	
		if(x$corr!=FALSE) prob1 <- apply(part1,1,pmnorm, mean=rep(0,2),varcov=covar1)	# y_1=1; y_2=0; outcome 3
		if(x$corr==FALSE) prob1 <- apply(part1,1,pmnorm, mean=rep(0,2),varcov=covar1.nc)	# y_1=1; y_2=0; outcome 3
		if(x$corr!=FALSE) prob2 <- apply(part2,1,pmnorm, mean=rep(0,2),varcov=covar2)	# y_1=1; y_2=1; outcome 4
		if(x$corr==FALSE) prob2 <- apply(part2,1,pmnorm, mean=rep(0,2),varcov=covar2.nc)	# y_1=1; y_2=1; outcome 4
	
		print(names(prob3))
		out.prob <- cbind(prob3, prob1, prob2)

####
		
		pred.out <- out.prob	
		if(prob==TRUE) colnames(pred.out) <- c("Prob of outcome 1", "Prob of outcome 3", "Prob of outcome 4")
		if(prob==FALSE) names(pred.out) <- c("Prob of outcome 1", "Prob of outcome 3", "Prob of outcome 4")
		
		yhat <- rep(NA,dim.x11[1])
		yhat[pred.out[,1]>0.5] <- 1
		yhat[pred.out[,1]<0.5& pred.out[,2]>pred.out[,3]] <- 3
		yhat[pred.out[,1]<0.5& pred.out[,3]>pred.out[,2]] <- 4  				
	 	}
	 	
	 	
	 	# code for a specific profile
	 if(missing(profile)==FALSE){
	  	dim.x11 <- x$DIM[c(1:2)]
 		dim.x14 <- x$DIM[c(3:4)]
 		dim.x24 <- x$DIM[c(5:6)]
 	
		X11 <- c(1,profile[1:(dim.x11[2]-1)])
		X14 <- c(1,profile[(dim.x11[2]):(dim.x14[2]-1)])
		X24 <- c(1,profile[(dim.x11[2]+dim.x14[2]-1):(dim.x11[2]+dim.x14[2]+dim.x24[2]-3)])
			
	 	beta <- x$coefficient
	 	# predicting outcomes
		a <- dim.x11[2]+ dim.x14[2] +1
		b <- length(beta)-1
		if (x$corr==FALSE) b <- length(beta)
		
		y24.latent <- t(as.matrix(X24))%*%as.matrix(beta[a:b])
		p24 <- pnorm(y24.latent)
		a <- dim.x11[2]+1
		b <- dim.x11[2]+dim.x14[2] 
		y14.latent <- p24*t(as.matrix(X14))%*%as.matrix(beta[a:b]) - t(as.matrix(X11))%*%as.matrix(beta[1:dim.x11[2]])
	
		rho <- beta[length(beta)]

###
		part1 <- cbind(-y24.latent,y14.latent); part2 <- cbind(y24.latent,y14.latent)		# combining for mvnorm
		covar1 <- matrix(c(1,-rho,-rho,1),2,2)				# negative rho
		covar2 <- matrix(c(1,rho,rho,1),2,2)	
		covar1.nc <- matrix(c(1,0,0,1),2,2)				# negative rho
		covar2.nc <- matrix(c(1,0,0,1),2,2)	
	
		prob3 <- 1-pnorm(y14.latent)	#  outcome 1	
		if(x$corr!=FALSE) prob1 <- apply(part1,1,pmnorm, mean=rep(0,2),varcov=covar1)	# y_1=1; y_2=0; outcome 3
		if(x$corr==FALSE) prob1 <- apply(part1,1,pmnorm, mean=rep(0,2),varcov=covar1.nc)	# y_1=1; y_2=0; outcome 3
		if(x$corr!=FALSE) prob2 <- apply(part2,1,pmnorm, mean=rep(0,2),varcov=covar2)	# y_1=1; y_2=1; outcome 4
		if(x$corr==FALSE) prob2 <- apply(part2,1,pmnorm, mean=rep(0,2),varcov=covar2.nc)	# y_1=1; y_2=1; outcome 4
	
		out.prob <- cbind(prob3, prob1, prob2)

####

#		rr <- beta[length(beta)]
#		vvv <- matrix(c(1,rr,rr,1),2,2)
#				p14 <- pmnorm(c(y14.latent,y24.latent),mean=c(0,0),varcov=vvv)	
	
#		pred.out <- rep(NA,3)
#		pred.out[1] <- 1-p14
#		pred.out[2] <- p14*(1-p24)
#		pred.out[3] <- p14*p24
	
		pred.out <- out.prob	
		if(prob==TRUE) colnames(pred.out) <- c("Prob of outcome 1", "Prob of outcome 3", "Prob of outcome 4")
		if(prob==FALSE) names(pred.out) <- c("Prob of outcome 1", "Prob of outcome 3", "Prob of outcome 4")
		
#		pXXX <- cbind(1-p14,p24)
#		colnames(pXXX) <- c("Prob of Player 1 going right", "Prob of Player 2 going right")
	
		yhat <- NA
		yhat[pred.out[1]>0.5] <- 1
		yhat[pred.out[1]<0.5& pred.out[2]>pred.out[3]] <- 3
		yhat[pred.out[1]<0.5& pred.out[2]<pred.out[3]] <- 4  					 	
	 	}
 	if(prob==FALSE) return(yhat)
 	if(prob==TRUE)  return(pred.out)
 }
