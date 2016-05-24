
boot_aldmck <- function(data, respondent=0, missing=NULL, polarity, iter=100){

	original <- aldmck(data=data, respondent=respondent, polarity=polarity,
		missing=missing,verbose=FALSE)
	samples <- matrix(NA, nrow=iter, ncol=length(original$stimuli))
	samples[1,] <- original$stimuli
	colnames(samples) <- names(original$stimuli)
	rownames(samples) <- 1:iter

	for(i in 2:iter){
        	tmpdat <- data[sample(1:nrow(data),nrow(data), replace=TRUE),]
		samples[i,] <- aldmck(data=tmpdat, respondent=respondent,
		polarity=polarity, missing=missing,verbose=FALSE)$stimuli
	}

	class(samples) <- "boot_aldmck"
	samples <- sign(samples[,1])*samples*-1	# For polarity consistency
	return(samples)
}

plot.boot_aldmck <- function(x, ...){

	if (class(x) != "boot_aldmck") stop("Data is not of class boot_aldmck.")

	x <- x[,order(colMeans(x))]
	positions <- colMeans(x)
	xrange <- ceiling(max(abs(positions))*10)/10
	plot(positions, 1:ncol(x), type="n", xlim=c(-xrange, xrange), yaxt="n",bty="n",
		ylim=c(0.8,ncol(x)),xlab="AM position score",ylab="",cex.lab=1.2, ...)
	bars <- apply(x,2,quantile,c(0.025,0.975))
	for(i in 1:ncol(bars)) segments(bars[1,i], i, bars[2,i],col="grey",lwd=2)
	points(positions, 1:ncol(bars), pch=20, cex=0.6)
	mtext(colnames(x), side=2, at=1:ncol(x),las=1,adj=1,cex=0.8, font=3)

}

plot.blackbox <- function(x, ...){

	if (class(x) != "blackbox") stop("Data is not of class blackbox.")

	hist(x$individuals[[x$dims]]$c1, breaks = seq(-1, 1, 0.1),
		main="Distribution of Blackbox Intercepts",
		xlab="Blackbox Intercept", ...)

}

plot.boot_blackbt <- function(x, ...){

	if (class(x) != "boot_blackbt") stop("Data is not of class boot_blackbt.")

	x <- x[,order(colMeans(x))]
	positions <- colMeans(x)
	xrange <- ceiling(max(abs(positions))*10)/10
	plot(positions, 1:ncol(x), type="n", xlim=c(-xrange, xrange), yaxt="n",bty="n",
		ylim=c(0.8,ncol(x)),xlab="Blackbox Transpose Score",ylab="",cex.lab=1.2, ...)
	bars <- apply(x,2,quantile,c(0.025,0.975))
	for(i in 1:ncol(bars)) segments(bars[1,i], i, bars[2,i],col="grey",lwd=2)
	points(positions, 1:ncol(bars), pch=20, cex=0.6)
	mtext(colnames(x), side=2, at=1:ncol(x),las=1,adj=1,cex=0.8, font=3)

}

boot_blackbt <- function(data, missing=NULL, dims=1, dim.extract=dims, minscale, iter=100){

	if(dim.extract > dims) stop("dim.extract must be less than dims")

	original <- blackbox_transpose(data=data, missing=missing, dims=dims, minscale=minscale, verbose=FALSE)
	samples <- matrix(NA, nrow=iter, ncol=nrow(original$stimuli[[dims]]))
	getdim <- paste("coord", dim.extract, "D",sep="")
	samples[1,] <- unlist(original$stimuli[[dims]][getdim])
	colnames(samples) <- rownames(original$stimuli[[dims]])
	rownames(samples) <- 1:iter

	for(i in 2:iter){
        	tmpdat <- data[sample(1:nrow(data),nrow(data), replace=TRUE),]
		rownames(tmpdat) <- NULL
		tmpres <- blackbox_transpose(tmpdat, missing=missing, dims=dims,
			minscale=minscale, verbose=FALSE)
		samples[i,] <- unlist(tmpres$stimuli[[dims]][getdim])

	        if(i %% 10 ==0){
			cat("\n\t\tIteration", i, "complete...")
			flush.console()
		}
	}

	class(samples) <- "boot_blackbt"
	samples <- sign(samples[,1])*samples*-1	# For polarity consistency
	return(samples)
}

predict.blackbox <- function(object, dims=1, ...){

	## Error catch for dims
    	if(!class(object)=="blackbox") stop("Input is not of class 'blackbox'.")
        if(dims < 1)  stop("dims must be great than 1")
	if(dims > object$dims)  stop(paste("dims must be equal or less than the number of estimate dimensions, which is", object$dims))

	## Extract object output
	W.hat<- as.matrix(object$stimuli[[dims]][,paste("w", 1:dims, sep="")])
	Psi.hat <- as.matrix(object$individuals[[dims]][,paste("c", 1:dims, sep="")])
	Jn <- rep(1, nrow(Psi.hat))
	c <- object$stimuli[[dims]][,"c"]

	## In blackbox, Keith already postmultiplied by sqrt(singular), so no need to do it here
	## This is NOT the case in blackbox_transpose, which needs to be multiplied by sqrt(singular)
	X.hat <- Psi.hat %*% t(W.hat) + Jn %o% c

	colnames(X.hat) <- rownames(object$stimuli[[dims]])
	rownames(X.hat) <- rownames(object$individuals[[dims]])
	return(X.hat)
}

predict.aldmck <- function(object, caliper=0.2, ...){

      	if(!class(object)=="aldmck") stop("Input is not of class 'aldmck'.")
      	if(caliper < 0) stop("Caliper must be positive.")
	N <-nrow(object$respondents)	#number of respondents
	j <- length(object$stimuli)	#number of stimuli
	Y <- matrix(rep(object$stimuli,N), nrow=N, ncol=j, byrow=TRUE)
	c <- matrix(object$respondent[,"intercept"],ncol=1) %*% matrix(rep(1,j),nrow=1)
	w <- matrix(object$respondent[,"weight"],ncol=1) %*% matrix(rep(1,j),nrow=1)
	X.hat <- (Y - c)/w
	dumpthese <- which(abs(object$respondent[,"weight"]) < caliper)
        X.hat[dumpthese,] <- NA

	colnames(X.hat) <- names(object$stimuli)
	rownames(X.hat) <- rownames(object$respondents)

	return(X.hat)
}

predict.blackbt <- function(object, dims=1, ...){

	## Error catch for dims
    	if(!class(object)=="blackbt") stop("Input is not of class 'blackbt'.")
	if(dims < 1)  stop("dims must be great than 1")
	if(dims > object$dims)  stop(paste("dims must be equal or less than the number of estimate dimensions, which is", object$dims))

	W.hat<- as.matrix(object$individuals[[dims]][,paste("w", 1:dims, sep="")])
	Psi.hat <- as.matrix(object$stimuli[[dims]][,paste("coord", 1:dims, "D", sep="")])
	Jn <- matrix(rep(1, nrow(Psi.hat)),ncol=1)
	c <- matrix(object$individuals[[dims]][,"c"],nrow=1)

	## In blackbox, Keith already postmultiplied by sqrt(singular), so no need to do it here
	## This is not the case in blackbox_transpose, which needs to be multiplied by sqrt(singular)
	singular <- sqrt(diag(object$fits$singular[1:dims]))
	X.hat <-  Psi.hat %*% singular %*% t(W.hat) + Jn %*% c
	X.hat <- t(X.hat)

	colnames(X.hat) <- rownames(object$stimuli[[dims]])
	rownames(X.hat) <- rownames(object$individuals[[dims]])

	return(X.hat)
}

summary.blackbox <- function(object, ...){

	x <- object
	if(!class(x)=="blackbox") stop("Input is not of class 'blackbox'.")
	cat("\n\nSUMMARY OF BLACKBOX OBJECT")
	cat("\n----------------------------------")
	for(i in 1:x$dims){
	  cat("\n")
	  print(x$stimuli[[i]], ...)
	}

	cat("\n\tDimensions Estimated:", x$dims)
	cat("\n\tNumber of Rows:", x$Nrow)
	cat("\n\tNumber of Columns:", x$Ncol)
	cat("\n\tTotal Number of Data Entries:", x$Ndata)
	cat("\n\tNumber of Missing Entries:", x$Nmiss)
	cat("\n\tPercent Missing Data: ", round(100*x$Nmiss/(x$Ndata + x$Nmiss),2),"%",sep="")
	cat("\n\tSum of Squares (Grand Mean):", x$SS_mean)
	cat("\n\n")
}

blackbox <- function(data,missing=NULL,verbose=FALSE,dims=1,minscale){

	### Error check each argument ###
	if(class(data) == "data.frame") data <- as.matrix(data)
	if(class(data) != "matrix") stop("Data is not a matrix, or data frame.")

	if(typeof(data) != "double") stop("Data are not numeric values, please convert it using as.numeric() or mode().")
	if(!is.null(missing) & !(is.matrix(missing) | is.vector(missing))) stop("Argument 'missing' must be a vector or matrix.")
	if(mode(missing) != "numeric" & !is.null(missing)) stop("Argument 'missing' must only contain integers.")
	if(!is.logical(verbose)) stop("Argument 'verbose' must be set TRUE or FALSE.")
	if(minscale<1) stop("Argument 'minscale' must be positive.")
	if(dims<1) stop("Argument 'dims' must be positive.")

	### Format the data input ###
	N <- nrow(data)
	NQ <- ncol(data)

	if(is.vector(missing))	data[data %in% missing] <- NA
	if(is.matrix(missing))	for(i in 1:ncol(data))	data[data[,i] %in% missing[,i],i] <- NA
	missval <- max(data,na.rm=TRUE) + 1	#code for missing data
	rawdata <- as.numeric(t(data))
	rawdata[is.na(rawdata)] <- missval

	##Longer output
	stimnames <- colnames(data)
	if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
	if(verbose){

	deleted <- sum(is.na(apply(data,1,sum)))
	cat("\n\n\tBeginning Blackbox Scaling...")

	#cat(nrow(data)-deleted, "of", nrow(data), "observations are complete.\n\t\t")
	cat(NQ, "stimuli have been provided.")

	}

	res <- .Fortran("blackbox",
                as.integer(N),			# NRESPONDENTS
                as.integer(NQ),			# NISSUES
		as.integer(dims),		# NDIMENSIONS, check later about mods
                as.integer(1),  		# NMISSING
	        as.numeric(rep(missval,NQ)),	# KMISS
		as.integer(minscale),		# MINSCALE
		as.integer(rep(1,N)),		# MID
                as.numeric(rawdata),		# KISSUE
		as.character("a"),		# CAND, not used
                fits = double(7*dims), 		# FITS
                psimatrix = double(N*((dims*(dims+1))/2)+2*N*dims),		# PSIMATRIX
                wmatrix = double((NQ)*((dims*(dims+1))/2)+2*(NQ)*dims),		# WMATRIX
		lresp = integer(N+NQ),						# LRESPONDENTS
		lmark = integer(N),		# LMARK
		fits2 = double(6),		# FITS2
                exitstatus = integer(1))	# EXITSTATUS

 	if (res$exitstatus != 1) stop("\n\n\t====== Blackbox did not execute properly ======\n\n")

	stimuli <- vector("list", dims)
	start <- 1
	end <- 3*NQ
	for(i in 1:dims){
		stimuli[[i]] <- as.data.frame(matrix(round(res$wmatrix[start:end],digits=3),nrow=NQ,ncol=i+2,byrow=T))
		colnames(stimuli[[i]]) <- c("c",paste("w",1:i,sep=""),"R2")
		rownames(stimuli[[i]]) <- stimnames
		stimuli[[i]] <- cbind(N=res$lresp[1:NQ],stimuli[[i]])
		start <- end + 1
		end <- start + (i+3)*NQ - 1
	}

	individuals <- vector("list", dims)
	start <- 1
	end <- N
	for(i in 1:dims){
		individuals[[i]] <- as.data.frame(matrix(round(res$psimatrix[start:end],digits=3),nrow=N,ncol=i,byrow=T))
		dumpthese <- (rowSums(individuals[[i]]==0)==i)
		individuals[[i]][dumpthese,] <- NA
		colnames(individuals[[i]]) <- c(paste("c",1:i,sep=""))
		rownames(individuals[[i]]) <- rownames(data)
		start <- end + 1
		end <- start + (i+1)*N - 1
	}

	fits <- matrix(res$fits,nrow=dims,ncol=7,byrow=T)
	fits <- as.data.frame(fits[,c(1:3,6:7),drop=FALSE])
	colnames(fits) <- c("SSE","SSE.explained","percent","SE","singular")
	rownames(fits) <- paste("Dimension",1:dims)

	result <- list(stimuli = stimuli, individuals=individuals, fits=fits,
		Nrow = res$fits2[1],
		Ncol = res$fits2[2],
		Ndata = res$fits2[3],
		Nmiss = res$fits2[4],
		SS_mean = res$fits2[6],
		dims=dims)

	class(result) <- c("blackbox")
	if (verbose) cat("\n\n\tBlackbox estimation completed successfully.\n\n")
	result

}

blackbox_transpose <- function(data,missing=NULL,verbose=FALSE,dims=1,minscale){

	### Error check each argument ###
	if(class(data) == "data.frame") data <- as.matrix(data)
	if(class(data) != "matrix") stop("Data is not a matrix, or data frame.")

	if(typeof(data) != "double") stop("Data are not numeric values, please convert it using as.numeric().")
	if(!is.null(missing) & !(is.matrix(missing) | is.vector(missing))) stop("Argument 'missing' must be a vector or matrix.")
	if(mode(missing) != "numeric" & !is.null(missing)) stop("Argument 'missing' must only contain integers.")
	if(!is.logical(verbose)) stop("Argument 'verbose' must be set TRUE or FALSE.")
	if(minscale<1) stop("Argument 'minscale' must be positive.")
	if(dims<1) stop("Argument 'dims' must be positive.")
	if(nrow(data)>1500) stop("There are more than N = 1500 respondents in the data.")

	### Format the data input ###
	N <- nrow(data)
	NQ <- ncol(data)

	if(is.vector(missing))	data[data %in% missing] <- NA
	if(is.matrix(missing))	for(i in 1:ncol(data))	data[data[,i] %in% missing[,i],i] <- NA
	missval <- max(data,na.rm=TRUE) + 1	#code for missing data
	rawdata <- as.numeric(t(data))
	rawdata[is.na(rawdata)] <- missval

	##Longer output
	stimnames <- colnames(data)
	if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
	if(verbose){

	deleted <- sum(is.na(apply(data,1,sum)))
	cat("\n\n\tBeginning Blackbox Transpose Scaling...")

	#cat(nrow(data)-deleted, "of", nrow(data), "observations are complete.\n\t\t")
	cat(NQ, "stimuli have been provided.")

	}

	res <- .Fortran("blackboxt",
                as.integer(N),			# NRESPONDENTS
                as.integer(NQ),			# NISSUES
		as.integer(dims),		# NDIMENSIONS, check later about mods
                as.integer(1),  		# NMISSING
	        as.double(rep(missval,NQ)),	# KMISS
		as.integer(minscale),		# MINSCALE
		as.integer(rep(1,N)),		# MID
                as.double(rawdata),		# KISSUE
		as.character("a"),		# CAND, not used
                fits = double(7*dims), 		# FITS
                psimatrix = double(N*((dims*(dims+1))/2)+2*N*dims),		# PSIMATRIX
                wmatrix = double((NQ)*((dims*(dims+1))/2)+2*(NQ)*dims),		# WMATRIX
		lresp = integer(N+NQ),						# LRESPONDENTS
		lmark = integer(N),		# LMARK
		fits2 = double(6),		# FITS2
                exitstatus = integer(1))	# EXITSTATUS

 	if (res$exitstatus != 1) stop("\n\n\t====== Blackbox-Transpose did not execute properly ======\n\n")

	stimuli <- vector("list", dims)
	start <- 1
	end <- 2*NQ
	for(i in 1:dims){
		stimuli[[i]] <- as.data.frame(matrix(round(res$wmatrix[start:end],digits=3),nrow=NQ,ncol=i+1,byrow=T))
		colnames(stimuli[[i]]) <- c(paste("coord",1:i,"D",sep=""),"R2")
		rownames(stimuli[[i]]) <- stimnames
		stimuli[[i]] <- cbind(N=res$lresp[(length(res$lresp)-NQ+1):length(res$lresp)],stimuli[[i]])
		start <- end + 1
		end <- start + (i+2)*NQ - 1
	}

	individuals <- vector("list", dims)
	start <- 1
	end <- 3*N
	for(i in 1:dims){
		individuals[[i]] <- as.data.frame(matrix(round(res$psimatrix[start:end],digits=3),nrow=N,ncol=i+2,byrow=T))
		colnames(individuals[[i]]) <- c("c",paste("w",1:i,sep=""),"R2")
		if(!is.null(rownames(data))) rownames(individuals[[i]]) <- rownames(data)
		start <- end + 1
		end <- start + (i+3)*N - 1
		individuals[[i]][!res$lmark,] <- NA
	}

	fits <- matrix(res$fits,nrow=dims,ncol=7,byrow=T)
	fits <- as.data.frame(fits[,c(1:3,6:7),drop=FALSE])
	colnames(fits) <- c("SSE","SSE.explained","percent","SE","singular")
	rownames(fits) <- paste("Dimension",1:dims)

	result <- list(stimuli = stimuli, individuals=individuals, fits=fits,
		Nrow = res$fits2[1],
		Ncol = res$fits2[2],
		Ndata = res$fits2[3],
		Nmiss = res$fits2[4],
		SS_mean = res$fits2[6],
		dims=dims)

	class(result) <- c("blackbt")
	if (verbose) cat("\n\n\tBlackbox-Transpose estimation completed successfully.\n\n")
	result
}

plot.blackbt <- function(x, xlim=c(-1,1), ...){

	if(!class(x)=="blackbt") stop("Input is not of class 'blackbt'.")

	colchoice <- rep(palette(),3)
	dens <- density(x$individuals[[1]]$w1,na.rm=TRUE)
	ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

	plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
		ylab="Density", main="Stimuli and Population Distribution", lab=c(5,5,7),
		col="blue", cex.main=1.2, cex.lab=1.2)

	for(i in 1:nrow(x$stimuli[[1]])){
		arrows(x$stimuli[[1]]$coord1D[i], 0.08, x$stimuli[[1]]$coord1D[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
		text(x=x$stimuli[[1]]$coord1D[i], y=0.1, rownames(x$stimuli[[1]])[i], srt=90, adj=0, font=2)
	}

	text(x=0.85*xlim[2],y=0.85*ymax,paste("N =",x$Ncol),cex=1.3)
}


plotcdf.blackbt <- function(x, align=NULL, xlim=c(-1.2,1), ...){

	if(!class(x)=="blackbt") stop("Input is not of class 'blackbt'.")

	colchoice <- rep(palette(),3)
	cdf <- ecdf(x$individuals[[1]]$w1)
	plot(cdf, lwd=2, xlim=xlim, lab=c(5,5,7), bty="n",
		cex.points=0.5, cex.main=1.2, cex.lab=1.2,
		xlab="Location", ylab="Cumulative Density",
		main="Stimuli and Population CDF", ...)

	abline(v=x$stimuli[[1]]$coord1D,lty=2,col='gray70')

	if(is.null(align)) align <- (min(cdf(x$stimuli[[1]]$coord1D))+xlim[1])/1.5
	for(i in 1:nrow(x$stimuli[[1]])){
		arrows(align, cdf(x$stimuli[[1]]$coord1D[i]), x$stimuli[[1]]$coord1D[i]-0.1, cdf(x$stimuli[[1]]$coord1D[i]), length=0.1, angle=20, col=colchoice[i], lwd=2)
	text(x=align, y=cdf(x$stimuli[[1]]$coord1D[i]), rownames(x$stimuli[[1]])[i], adj=1, font=2)
	}
}

summary.blackbt <- function(object, ...){

	x <- object
	if(!class(x)=="blackbt") stop("Input is not of class 'blackbt'.")
	cat("\n\nSUMMARY OF BLACKBOX TRANSPOSE OBJECT")
	cat("\n----------------------------------")
	for(i in 1:x$dims){
	  cat("\n")
	  print(x$stimuli[[i]],...)
	}

	cat("\n\tDimensions Estimated:", x$dims)
	cat("\n\tNumber of Rows:", x$Nrow)
	cat("\n\tNumber of Columns:", x$Ncol)
	cat("\n\tTotal Number of Data Entries:", x$Ndata)
	cat("\n\tNumber of Missing Entries:", x$Nmiss)
	cat("\n\tPercent Missing Data: ", round(100*x$Nmiss/(x$Ndata + x$Nmiss),2),"%",sep="")
	cat("\n\tSum of Squares (Grand Mean):", x$SS_mean)
	cat("\n\n")
}

plot.aldmck <- function(x, ...){

	if(!class(x)=="aldmck") stop("Input is not of class 'aldmck'.")

        op <- par(no.readonly=TRUE)     
	par(mfrow=c(2,2))
	plot.AM(x,...)
	plot.cdf(x,...)
	plot.aldmck_positive(x,...)
	plot.aldmck_negative(x,...)
        suppressWarnings(par(op))

}

plot.AM <- function(x, xlim=c(-2,2), ...){

	if(!class(x)=="aldmck") stop("Input is not of class 'aldmck'.")

	if(sum(x$respondents[,"weight"]>0,na.rm=TRUE)==0){
		plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
		text(0.5,0.5,"No\nself\nplacement",cex=3)
		return()
	}

	colchoice <- rep(palette(),3)
	dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
	ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

	plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
		ylab="Density", main="Stimuli and Population Distribution", lab=c(5,5,7),
		col="blue", cex.main=1.2, cex.lab=1.2, ...)

	for(i in 1:length(x$stimuli)){
		arrows(x$stimuli[i], 0.08, x$stimuli[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
		text(x=x$stimuli[i], y=0.1, names(x$stimuli)[i], srt=90, adj=0, font=2)
	}

	text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N),cex=1.3)
}


plot.aldmck_positive <- function(x, xlim=c(-2,2), ...){

	if(!class(x)=="aldmck") stop("Input is not of class 'aldmck'.")

	if(sum(x$respondents[,"weight"]>0,na.rm=TRUE)==0){
		plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
		text(0.5,0.5,"No\nweights",cex=3)
		return()
	}

	x$respondents <- x$respondents[x$respondents[,"weight"]>0,]

	colchoice <- rep(palette(),3)
	dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
	ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

	plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
		ylab="Density", main="Positive Weights", lab=c(5,5,7),
		col="blue", cex.main=1.2, cex.lab=1.2, ...)

	for(i in 1:length(x$stimuli)){
		arrows(x$stimuli[i], 0.08, x$stimuli[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
		text(x=x$stimuli[i], y=0.1, names(x$stimuli)[i], srt=90, adj=0, font=2)
	}

	text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N.pos),cex=1.3)
}

plot.aldmck_negative <- function(x, xlim=c(-2,2), ...){

	if(!class(x)=="aldmck") stop("Input is not of class 'aldmck'.")

	if(sum(x$respondents[,"weight"]<0,na.rm=TRUE)==0){
		plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
		text(0.5,0.5,"No\nnegative\nweights",cex=2.5)
		return()
	}

	x$respondents <- x$respondents[x$respondents[,"weight"]<0,]
	dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
	ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

	if(x$N.neg>50){
	colchoice <- rep(palette(),3)
	dens <- density(x$respondents[,"idealpt"],na.rm=TRUE)
	ymax <- max(dens$y) - max(dens$y) %% 0.2 + 0.2

	plot(dens, xlim=xlim, ylim=c(0,ymax), lwd=3, bty="n", xlab="Location",
		ylab="Density", main="Negative Weights", lab=c(5,5,7),
		col="blue", cex.main=1.2, cex.lab=1.2, ...)

	for(i in 1:length(x$stimuli)){
		arrows(x$stimuli[i], 0.08, x$stimuli[i], 0, length=0.1, angle=20, col=colchoice[i], lwd=2)
		text(x=x$stimuli[i], y=0.1, names(x$stimuli)[i], srt=90, adj=0, font=2)
	}
	text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N.neg),cex=1.3)
	} #end if N>50

	if(x$N.neg<51){
	hist(x$respondents[,"idealpt"],xlim=c(-2,2),ylim=c(0,ymax),breaks=30,col="lightgrey",freq=F,
		xlab="Location",ylab="Density", main="Stimuli and Population Distribution: Negative Weights", cex.main=1.2, cex.lab=1.2,...)
	text(x=0.9*xlim[2],y=0.9*ymax,paste("N =",x$N.neg),cex=1.3)

	} #end if N<51
}

plot.cdf <- function(x, align=NULL, xlim=c(-2,2), ...){

	if(!class(x)=="aldmck") stop("Input is not of class 'aldmck'.")

	if(sum(x$respondents[,"weight"]>0,na.rm=TRUE)==0){
		plot(0,1,xlim=c(0,1),ylim=c(0,1),type="n",xaxt="n",yaxt="n",xlab="", ylab="",bty="n")
		text(0.5,0.5,"No\nself\nplacement",cex=3)
		return()
	}

	colchoice <- rep(palette(),3)
	cdf <- ecdf(x$respondents[,"idealpt"])
	plot(cdf, lwd=2, xlim=xlim, lab=c(5,5,7), bty="n",
		cex.points=0.5, cex.main=1.2, cex.lab=1.2,
		xlab="Location", ylab="Cumulative Density",
		main="Stimuli and Population CDF", ...)

	abline(v=x$stimuli,lty=2,col='gray70')

	if(is.null(align)) align <- cdf(x$stimuli[1]) - 1.5
	for(i in 1:length(x$stimuli)){
		arrows(align + 0.1, cdf(x$stimuli[i]), x$stimuli[i]-0.1, cdf(x$stimuli[i]), length=0.1, angle=20, col=colchoice[i], lwd=2)
		text(x=align, y=cdf(x$stimuli[i]), names(x$stimuli)[i], adj=1, font=2)
	}
}


aldmck <- function(data, respondent = 0, missing=NULL, polarity, verbose=FALSE) {

	### Error check each argument ###
	if(class(data) == "data.frame") data <- as.matrix(data)
	if(class(data) != "matrix") stop("Data is not a matrix or data frame.")
	if(typeof(data) != "double") stop("Data are not numeric values, please convert it using as.numeric() or mode().")

	if(mode(respondent) != "numeric")  stop("Respondent is not specified as an integer.")
	if(respondent > ncol(data))  stop("Respondent set to a column greater than number of columns in data.")
	if(respondent < 0)  stop("Respondent cannot be negative, set respondent=0 is self identification is unavailable.")

	if(mode(polarity) != "numeric")  stop("Polarity is not specified as an integer.")
	if(polarity == respondent) stop("Self-placements cannot also be set as polarity.")
	if(polarity > ncol(data))  stop("Polarity set to a column greater than number of columns in data.")
	if(polarity < 0)  stop("Polarity cannot be negative.")

	if(!is.null(missing) & !(is.matrix(missing) | is.vector(missing))) stop("Argument 'missing' must be a vector or matrix.")
	if(mode(missing) != "numeric" & !is.null(missing)) stop("Argument 'missing' must only contain integers.")

	if(!is.logical(verbose)) stop("Argument 'verbose' must be set TRUE or FALSE.")

	### Format the data input ###
	N <- nrow(data)
	NQ <- ncol(data) - 1
	NRESP <- respondent
	if(respondent==0){
		NQ <- ncol(data)
		NRESP  <- ncol(data) + 1
	}

	if(is.vector(missing))	data[data %in% missing] <- NA
	if(is.matrix(missing))	for(i in 1:ncol(data))	data[data[,i] %in% missing[,i],i] <- NA
	missval <- max(data,na.rm=TRUE) + 1	#code for missing data
	if(respondent==0) data <- cbind(data, fakedata = rep(missval,nrow(data)))
	rawdata <- as.integer(t(data))
	rawdata[is.na(rawdata)] <- missval

	##Longer output
	deleted <- sum(is.na(apply(data,1,sum)))
	stimnames <- colnames(data)
	if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
	if(verbose){

	cat("\n\n\tBeginning Aldrich-McKelvey Scaling...")
	if(respondent!=0) cat("\n\n\t\tColumn '",stimnames[respondent],"' is set as the self placement.", sep="")
	if(respondent==0) cat("\n\n\t\tSelf-placements have not been provided by the user.")

	cat("\n\t\tColumn '",stimnames[polarity],"' is set as the left-leaning stimulus.\n\t\t", sep="")
	cat(nrow(data)-deleted, "of", nrow(data), "observations are complete.\n\t\t")
	cat(NQ, "stimuli have been provided.")
	}
    if ((respondent != 0) & (respondent < polarity))   polarity = polarity - 1


	### Send the data to Fortran ###
	res <- .Fortran("mckalnew",
              as.integer(N), 			# NRESPONDENTS
              as.integer(NQ+1),	 		# NISSUES 
              as.integer(NRESP),  		# NSELFPOS
              as.integer(1),			# NMISSING (maximum number of missing valaues)
	      as.numeric(rep(missval,NQ+1)),	# KMISS, input number of missing values in each stimuli, length(NMISS*(NQ+1))), need to be double later
              as.integer(polarity),       	# POLARITY (input, which is on left)
              as.integer(rep(1,N)),		# MID, input, ID number of individual,length(nrow(data))   
              as.numeric(rawdata),		# KISSUE, data matrix, needs to be double later, length(nrow(data)*(NQ+1)))
#              as.character("a"),		# CAND, input, names of stimuli each column,length(NQ+1)
              fits = double(5),			# FITS (output), AMfit, R^2, numberscaled, negweights, posweights 
              psimatrix = double(N*4),		# PSIMATRIX, intercept C, W weight, scaleposition, R^2 
              stimcoord = double(NQ),		# STIMCOORDS, stimulus coordinates
              eigenvalues = double(NQ),		# EIGENVALUES
              exitstatus = integer(1))		# EXITSTATUS
	if(res$exitstatus!=1) stop("\n\n\t====== Aldrich-McKelvey did not execute properly ======\n\n")

	### Format the output ###
	stimnames <- colnames(data)
	if(is.null(stimnames)) stimnames <- paste("stim", 1:N, sep="")
	stimuli <- as.numeric(res$stimcoord)
	names(stimuli) <- stimnames[-NRESP]

	respondents <- matrix(res$psimatrix, ncol=4, byrow=T)
	respondents[respondents==0] <- NA
	midnames <- as.character(unlist(rownames(data)))
	if(is.null(midnames)) midnames <- paste("resp", 1:N, sep="")
	rownames(respondents) <- midnames
	colnames(respondents) <- c("intercept", "weight", "idealpt", "R2")

	## For JSS edit, R2 removed
	respondents <- respondents[,1:3]

    if (respondent == 0) {
        selfplace <- rep(NA, N)
	data <- data[,-ncol(data)]
    }    else {
        selfplace <- data[, respondent]
        data <- data[, -respondent]
    }
    
	polinfo <- suppressWarnings(apply(data,1,cor,y=stimuli,use="pairwise.complete.obs"))
	polinfo[which(respondents[,"weight"] < 0)] <- 0
	respondents <- cbind(respondents, selfplace, polinfo)


	result <- list(stimuli = stimuli, respondents = as.data.frame(respondents),
			eigenvalues = as.numeric(res$eigenvalues),
			AMfit = as.numeric(res$fits[1]), 
#			R2 = as.numeric(res$fits[2]),
#			N = nrow(data)-deleted, N.neg = as.integer(res$fits[5]),
			N = as.integer(res$fits[4])+as.integer(res$fits[5]), N.neg = as.integer(res$fits[5]),
			N.pos = as.integer(res$fits[4]))

	class(result) <- c("aldmck")
	if(verbose) cat("\n\n\tAldrich-McKelvey estimation completed successfully.\n\n")    
	result
 

}

summary.aldmck <- function(object, ...){

    x<-object
    if(!class(x)=="aldmck") stop("Input is not of class 'aldmck'.")

    cat("\n\nSUMMARY OF ALDRICH-MCKELVEY OBJECT")
    cat("\n----------------------------------")
    cat("\n\nNumber of Stimuli:", length(x$stimuli))
    cat("\nNumber of Respondents Scaled:", x$N)
     if(sum(!is.na(x$respondents[,"selfplace"])) != 0){
	cat("\nNumber of Respondents (Positive Weights):", x$N.pos)
	cat("\nNumber of Respondents (Negative Weights):", x$N.neg)
#	cat("\n\nR-Squared:", round(x$R2, digits=2)
     }
    cat("\nReduction of normalized variance of perceptions:", round(x$AMfit, digits=2), "\n\n")


    final <- matrix(round(sort(x$stimuli),3),ncol=1)
    rownames(final) <- names(sort(x$stimuli))
    colnames(final) <- c("Location")
    print(final,...)
    cat("\n\n")
}

stimuli <- function(object){
        if(class(object)=="blackbox" | class(object)=="blackbt") return(object$stimuli)
        if(class(object)=="aldmck") return(object$stimuli)
        stop("Input is not of class 'blackbox' or 'blackbt' or 'aldmck'.")
}

individuals <- function(object){
        if(class(object)=="blackbox" | class(object)=="blackbt") return(object$individuals)
        if(class(object)=="aldmck") return(object$respondents)
        stop("Input is not of class 'blackbox' or 'blackbt' or 'aldmck'.")
}

fit <- function(object){
        if(class(object)=="blackbox" | class(object)=="blackbt") return(object$fits)
        if(class(object)=="aldmck") return(object$AMfit)
        stop("Input is not of class 'blackbox' or 'blackbt' or 'aldmck'.")
}

dim.blackbox <- function(x){
	return(c(x$Nrow, x$Ncol))
}

dim.blackbt <- function(x){
	return(c(x$Nrow, x$Ncol))
}

dim.aldmck <- function(x){
	return(c(x$N, length(x$stimuli)))
}

ncol.blackbox <- function(x){
	return(x$Ncol)
}

nrow.blackbox <- function(x){
	return(x$Nrow)
}

ncol.blackbt <- function(x){
	return(x$Ncol)
}

nrow.blackbt <- function(x){
	return(x$Nrow)
}

nrow.aldmck <- function(x){
	return(x$N)
}

ncol.aldmck <- function(x){
	return(length(x$stimuli))
}




