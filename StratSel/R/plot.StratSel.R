plot.StratSel <-
function(x, profile, x.move, x.range, uncertainty=FALSE, n.sim=100, ci=0.95, ylim, xlab, ylab1, ylab2, ylab3, plot.nr, ...){
	object <- x
	if(missing(profile)) stop(gettextf("You have to specify a profile -- see help(plot.StratSel)."))
	if(missing(x.move)) stop(gettextf("You have to specify which X should vary in the plots -- see help(plot.StratSel)."))
	if(missing(x.range)) stop(gettextf("You have to specify the range of the moving X -- see help(plot.StratSel)."))
	test.length <- length(profile)
	#if(object$corr==FALSE) stop(gettextf("You are attempting to plot a model with uncorrelated errors. Use functionality of games package for this."))

	if(missing(ylim)) YLIM <- c(0,1)
	
	if(missing(xlab)) XLAB <- paste("X from", x.range[1], "to", x.range[2], seq=" ")
	if(missing(ylab1)) YLAB1 <- "Prob of Outcome 1"
	if(missing(ylab2)) YLAB2 <- "Prob of NOT Outcome 1"
	if(missing(ylab3)) YLAB3 <- "Prob of Outcome 4"
	if(missing(plot.nr)) plot.nr <- c(1,2,3)

	if(object$corr==FALSE){
		covar <- object$vcov		
		beta <- object$coefficient
		test.size <- test.length - length(beta)
		if(test.size!=0) stop(gettextf("You are not supplying a valid profile. Too many or too few values."))
		}

	if(object$corr!=FALSE){
		covar <- object$vcov[-dim(object$vcov)[1],-dim(object$vcov)[2]]		
		beta <- object$coefficient[-length(object$coefficient)]
		test.size <- test.length - length(beta)
		if(test.size!=0) stop(gettextf("You are not supplying a valid profile. Too many or too few values."))
		}

	if(uncertainty==TRUE){
		beta <- mvrnorm(n.sim, beta, covar)
		DIM <- object$DIM[c(2,4,6)]	 
		beta11 <- beta[,c(1:DIM[1])]
		beta14 <- beta[,c((DIM[1]+1):(DIM[1]+DIM[2]))]
		beta24 <- beta[,c((DIM[1]+DIM[2]+1):(DIM[1]+DIM[2]+DIM[3]))]
		DAT <- matrix(profile,n.sim,length(profile), byrow=TRUE)
		XX <- seq(from=x.range[1],to=x.range[2],length.out=n.sim)
		DAT[,x.move] <- XX
		x11 <- DAT[,c(1:DIM[1])]
		x14 <- DAT[,c((DIM[1]+1):(DIM[1]+DIM[2]))]
		x24 <- DAT[,c((DIM[1]+DIM[2]+1):(DIM[1]+DIM[2]+DIM[3]))]
		y24.latent <- x24%*%t(beta24)
		P24 <- t(apply(y24.latent,1,pnorm))
		y14.latent <- P24*x14%*%t(beta14) - x11%*%t(beta11)
		P14 <- t(apply(y14.latent,1,pnorm))
		p14.m <- rowMeans(P14)
		p24.m <- rowMeans(P24)
		p14.sort <- t(apply(P14,1,sort))
		p24.sort <- t(apply(P24,1,sort))
		low <- (1-ci)/2
		high <- 1- (1-ci)/2
		L <- round(n.sim*low); 	H <- round(n.sim*high)
		p14.l <- p14.sort[,L]
		p14.h <- p14.sort[,H]
		p24.l <- p24.sort[,L]
		p24.h <- p24.sort[,H]		

		if (1 %in% plot.nr)    plot(y=1-p14.m, XX, ylab=YLAB1, type="l", xlab=XLAB, ylim=YLIM, ...)
		if (1 %in% plot.nr) 	points(XX, 1-p14.l, type="l", lty=2)	
		if (1 %in% plot.nr) 	points(XX, 1-p14.h, type="l", lty=2)	
		if (2 %in% plot.nr) 	plot(y=p14.m, XX, ylab=YLAB2, type="l", ylim=YLIM, xlab=XLAB, ...)	
		if (2 %in% plot.nr) 	points(XX, p14.l, type="l", lty=2)	
		if (2 %in% plot.nr) 	points(XX, p14.h, type="l", lty=2)	
		if (3 %in% plot.nr) 	plot(y=p24.m, XX, ylab=YLAB3, type="l", ylim=YLIM, xlab=XLAB, ...)	
		if (3 %in% plot.nr) 	points(XX, p24.l, type="l", lty=2)	
		if (3 %in% plot.nr) 	points(XX, p24.h, type="l", lty=2)	
		} 

	if(uncertainty==FALSE){
		DIM <- object$DIM[c(2,4,6)]	 
		beta11 <- beta[c(1:DIM[1])]
		beta14 <- beta[c((DIM[1]+1):(DIM[1]+DIM[2]))]
		beta24 <- beta[c((DIM[1]+DIM[2]+1):(DIM[1]+DIM[2]+DIM[3]))]
		DAT <- matrix(profile,n.sim,length(profile), byrow=TRUE)
		XX <- seq(from=x.range[1],to=x.range[2],length.out=n.sim)
		DAT[,x.move] <- XX
		x11 <- DAT[,c(1:DIM[1])]
		x14 <- DAT[,c((DIM[1]+1):(DIM[1]+DIM[2]))]
		x24 <- DAT[,c((DIM[1]+DIM[2]+1):(DIM[1]+DIM[2]+DIM[3]))]
		y24.latent <- x24%*%(as.matrix(beta24))
		P24 <- t(apply(y24.latent,1,pnorm))
		y14.latent <- x14%*%(beta14) * t(P24) - x11%*%(beta11)
		P14 <- t(apply(y14.latent,1,pnorm))

		if (1 %in% plot.nr) plot(y=1-P14, XX, ylab=YLAB1, type="l", xlab=XLAB, ylim=YLIM, ...)
		if (2 %in% plot.nr) plot(y=P14, XX, ylab=YLAB2, type="l", ylim=YLIM, xlab=XLAB, ...)	
		if (3 %in% plot.nr) plot(y=P24, XX, ylab=YLAB3, type="l", ylim=YLIM, xlab=XLAB, ...)	
		} 
}
