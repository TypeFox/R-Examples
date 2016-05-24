LSCV.density <- function(data, hlim = NULL, res = 128, edge = TRUE, WIN = NULL, quick = TRUE, comment = TRUE){
	if(comment) cat("Initialising...\n")
	if(!is.null(hlim)) if(hlim[1] >= hlim[2]) stop("invalid h limits")
	if(class(data)=="ppp"){
		data.ppp <- data
	} else if(class(data)=="list"){
		if(is.null(WIN)) stop("must supply WIN if data is not a ppp.object")
		if(is.null(data$x)||is.null(data$y)) stop("data list must have two components named 'x' and 'y'")
		if(length(data$x)!=length(data$y)) stop("data components 'x' and 'y' are of unequal lengths")
		if(length(data$x)<10) warning("less than 10 case observations!")
		data.ppp <- ppp(x=data$x,y=data$y,window=WIN)
	} else if((class(data)=="matrix")||(class(data)=="data.frame")){
		if(is.null(WIN)) stop("must supply WIN if data is not a ppp.object")
		if(length(na.omit(data[,1]))!=length(na.omit(data[,2]))) stop("data components 'x' and 'y' are of unequal lengths in 'data'")
		if(nrow(data)<10) warning("less than 10 case observations!")
		data.ppp <- ppp(x=data[,1],y=data[,2],window=WIN)
	} else {
		stop("invalid data object")
	}
	
	W <- data.ppp$window
	M <- N <- res+1
	m <- 0:(M-1)
	n <- 0:(N-1)
	del1 <- (W$xrange[2]-W$xrange[1])/(M-1)
	del2 <- (W$yrange[2]-W$yrange[1])/(N-1) 
	mgr <- W$xrange[1]+m*del1
	ngr <- W$yrange[1]+n*del2
	mcens <- W$xrange[1]+.5*del1+(0:(M-2))*del1
	ncens <- W$yrange[1]+.5*del2+(0:(N-2))*del2
	cellinside.vec <- rep(0,(M-1)*(N-1))
	cellCentroids <- list(x=rep(NA,(M-1)*(N-1)),y=rep(NA,(M-1)*(N-1)))
	pos <- 0
	for(i in 1:(M-1)){     
		for(j in 1:(N-1)){
			temp.win <- owin(xrange=c(mgr[i],mgr[i+1]),yrange=c(ngr[j],ngr[j+1]))
			pos <- pos+1
			temp.cen <- centroid.owin(w=temp.win)
			cellCentroids$x[pos] <- temp.cen$x
			cellCentroids$y[pos] <- temp.cen$y
			if(inside.owin(x=temp.cen$x,y=temp.cen$y,w=W)) cellinside.vec[pos] <- 1
		}
	}
	cellinside.mat <- t(matrix(cellinside.vec,N-1,M-1))
	selector <- matrix(FALSE,M-1,N-1)
	selector[matrix(as.logical(cellinside.mat),M-1,N-1)] <- TRUE
	
	if(!is.null(data.ppp$marks)) data.ppp$marks <- NULL
	if(is.null(hlim)){
		md <- min(nndist(unique(data.ppp)))
		hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
	}
	
	if(comment) cat(paste("...searching for optimal h within (",hlim[1],",",hlim[2],")...\n",sep=""))
	if(quick){
		hopt <- optimise(LSCV.density.single,interval=hlim,data.ppp=data.ppp,res=res,edge=edge,inside=selector)$minimum
		if(comment) cat("...done.\n")
		return(hopt)
	} else {
		hseq <- seq(hlim[1],hlim[2],length=100)
		lscv.vec <- c()
		if(comment) pb <- txtProgressBar(min=0,max=length(hseq),style=3)
		for(i in 1:length(hseq)){
			lscv.vec[i] <- LSCV.density.single(hseq[i],data.ppp=data.ppp,res=res,edge=edge,inside=selector)
			if(comment) setTxtProgressBar(pb,i)
		}
		if(comment) close(pb)
		
		ind <- which(lscv.vec==min(lscv.vec[!is.na(lscv.vec)]))
		if((ind==1)||(ind==length(hseq))) warning("bandwidth limit selected -- try a different hlim")
		
		return(list(hopt=hseq[!is.na(lscv.vec)][which(lscv.vec[!is.na(lscv.vec)]==min(lscv.vec[!is.na(lscv.vec)]))],lscv=lscv.vec,ind=ind))
	}
}
