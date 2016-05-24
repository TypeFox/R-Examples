np.graph <-
function(object, which, var, exp, simul, obs, xlab, ylab, xlim, ylim, main){

	if(missingArg(simul)) simul <- TRUE
	if(missingArg(obs)) obs <- FALSE
	if(missingArg(exp)) exp <- FALSE
	if(missingArg(main) || !is.character(main)) main <- " "

	if(missingArg(which))
	stop("The value of the which argument is missing!!",call.=FALSE)

	if(which!=1 & which!=2)
	stop("value of the which argument is invalid!!",call.=FALSE)

	reem <- function(aa,b){
		 ag <- aa
		 ag <- sub("(", "", ag,fixed=TRUE)
		 ag <- sub(")", "", ag,fixed=TRUE)
		 ag <- sub(b, "", ag,fixed=TRUE)
		 ag <- strsplit(ag, ",")
		 ag <- ag[[1]][1]
		 ag
	}

	if(which==1){
		if(sum(object$qm) == 0) stop("There are not nonparametric components in the Median/Location submodel!!",call.=FALSE)
	    xb <- colnames(object$model.matrix.mu)[(object$p+1):(length(object$qm)+object$p)]
		if(missingArg(var) && length(xb)==1) var <- xb
	    idx2 <- grepl(var,xb,fixed=TRUE)
		if(sum(idx2) == 0) stop("There is not the nonparametric effect requested by the user!!",call.=FALSE)
		varx2 <- min(seq(1,length(xb),by=1)[idx2==1])
		which.mu <- object$which.mu[object$which.mu != ""]
		posi <- c(0,cumsum(object$qm)) + object$p
		xx <- object$model.matrix.mu[,(object$p + varx2)]
		xxm <- eval(parse(text=sub(reem(xb[varx2],which.mu[varx2]),"xx",xb[varx2],fixed=TRUE)))
		cent <- qr.Q(qr(matrix(1,ncol(attr(xxm,"N")),1)),complete=TRUE)[,-1]
		psps <- attr(xxm,"N")%*%cent%*%object$theta.mu[(posi[varx2] + 1):posi[varx2+1]]
		psps2 <- attr(xxm,"N2")%*%cent%*%object$theta.mu[(posi[varx2] + 1):posi[varx2+1]]
		response <- object$z_es*sqrt(object$phi.fitted) + psps
		se <- sqrt(diag(attr(xxm,"N2")%*%cent%*%tcrossprod(object$vcov.mu[(posi[varx2] + 1):posi[varx2+1],(posi[varx2] + 1):posi[varx2+1]],attr(xxm,"N2")%*%cent)))
		if(!simul) nm <- 1 else nm <- length(as.numeric(levels(factor(xx))))
		if(missingArg(xlab) || !is.character(xlab))xlab <- reem(xb[varx2],which.mu[varx2])
 	    if(missingArg(ylab) || !is.character(ylab))ylab <- ifelse(exp,paste("exp(",xb[varx2],")"),xb[varx2])
	}
	if(which==2){
		if(sum(object$q) == 0) stop("There are not nonparametric components in the Skewness/Dispersion submodel!!",call.=FALSE)
	    xb <- colnames(object$model.matrix.phi)[(object$l+1):(length(object$q)+object$l)]
		if(missingArg(var) && length(xb)==1) var <- xb
	    idx2 <- grepl(var,xb,fixed=TRUE)
		if(sum(idx2) == 0) stop("There is not the nonparametric effect requested by the user!!",call.=FALSE)
		varx2 <- min(seq(1,length(xb),by=1)[idx2==1])
		which.phi <- object$which.phi[object$which.phi != ""]
		posi <- c(0,cumsum(object$q)) + object$l
		xx <- object$model.matrix.phi[,(object$l + varx2)]
		xxm <- eval(parse(text=sub(reem(xb[varx2],which.phi[varx2]),"xx",xb[varx2],fixed=TRUE)))
		cent <- qr.Q(qr(matrix(1,ncol(attr(xxm,"N")),1)),complete=TRUE)[,-1]
		psps <- attr(xxm,"N")%*%cent%*%object$theta.phi[(posi[varx2] + 1):posi[varx2+1]]
		response <- log((object$z_es*sqrt(object$phi.fitted))^2/object$xix) - log(object$phi.fitted) + psps
		psps2 <- attr(xxm,"N2")%*%cent%*%object$theta.phi[(posi[varx2] + 1):posi[varx2+1]]
		se <- sqrt(diag(attr(xxm,"N2")%*%cent%*%tcrossprod(object$vcov.phi[(posi[varx2] + 1):posi[varx2+1],(posi[varx2] + 1):posi[varx2+1]],attr(xxm,"N2")%*%cent)))
		if(!simul) nm <- 1 else nm <- length(as.numeric(levels(factor(xx))))
		if(missingArg(xlab) || !is.character(xlab))xlab <- reem(xb[varx2],which.phi[varx2])
 	    if(missingArg(ylab) || !is.character(ylab))ylab <- ifelse(exp,paste("exp(",xb[varx2],")"),xb[varx2])
	}

	  xx2 <- seq(min(xx),max(xx),length=200)
	  if(missingArg(xlim)) xlim <- range(xx)
	  if(exp==TRUE){
	  	  lims <- cbind(exp(psps2 - qnorm(0.025/nm)*se),exp(psps2 + qnorm(0.025/nm)*se))
	  	  if(missingArg(ylim)) if(obs) ylim <- exp(range(response)) else ylim <- range(lims)
		  plot(xx2, exp(psps2), xlim=xlim, ylim=ylim, type="l", col="white", xlab=xlab, ylab=ylab)
		  polygon(c(xx2,xx2[200:1]), c(lims[,1],lims[200:1,2]), col="light gray",border="light gray")
		  if(obs){
			  par(new=TRUE)
			  if(object$censored==FALSE)
    	      	plot(xx, exp(response), xlim=xlim, ylim=ylim, type="p", cex=0.3, lwd=3, xlab="", ylab="", main="")
			  else{
    	      	plot(xx[object$event==0], exp(response[object$event==0]), xlim=xlim, ylim=ylim, type="p", cex=0.3, lwd=3, xlab="", ylab="", main="")
				par(new=TRUE)
    	      	plot(xx[object$event==1], exp(response[object$event==1]), xlim=xlim, ylim=ylim, pch="+", xlab="", ylab="", main="")
			  }
		  }
		  par(new=TRUE)
		  plot(xx2, exp(psps2), xlim=xlim, ylim=ylim, type="l", col="blue", xlab=xlab, ylab=ylab, main=main)

	  }else{
	  	  lims <- cbind(psps2 - qnorm(0.025/nm)*se,psps2 + qnorm(0.025/nm)*se)
	  	  if(missingArg(ylim)) if(obs) ylim <- range(response) else ylim <- range(lims)
		  plot(xx2, psps2, xlim=xlim, ylim=ylim, type="l", col="white", xlab=xlab, ylab=ylab)
		  polygon(c(xx2,xx2[200:1]), c(lims[,1],lims[200:1,2]), col="light gray",border="light gray")
		  if(obs){
			  par(new=TRUE)
			  if(object$censored==FALSE)
		          plot(xx, response, xlim=xlim, ylim=ylim, type="p", cex=0.3, lwd=3, xlab="", ylab="", main="")
			  else{
		          plot(xx[object$event==0], response[object$event==0], xlim=xlim, ylim=ylim, type="p", cex=0.3, lwd=3, xlab="", ylab="", main="")
				  par(new=TRUE)
		          plot(xx[object$event==1], response[object$event==1], xlim=xlim, ylim=ylim, pch="+", xlab="", ylab="", main="")
			  }
		  }
		  par(new=TRUE)
		  plot(xx2, psps2, xlim=xlim, ylim=ylim, type="l", col="blue", xlab=xlab, ylab=ylab, main=main)
	  }
#	  list(xx2=xx2, psps2=psps2, lims=lims, xx=xx, response=response)
}
