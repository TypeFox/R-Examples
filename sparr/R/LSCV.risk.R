LSCV.risk <- function(cases, controls, hlim = NULL, method = c("kelsall-diggle", "hazelton"), res = 128, WIN = NULL, edge = TRUE, comment = TRUE){
	if(comment) cat("Initialising...\n")
	if(!is.null(hlim)) if(hlim[1] >= hlim[2]) stop("invalid h limits")
	if(class(cases)!=class(controls)) stop("case and control arguments must be of the same class")
	if((method[1]!="kelsall-diggle")&&(method[1]!="hazelton")) stop("method of selection must be either 'kelsall-diggle' or 'hazelton'")
	
	if(class(cases)=="ppp"){
		cases.ppp <- cases
		controls.ppp <- controls
	} else if(class(cases)=="list"){
		if(is.null(WIN)) stop("must supply WIN if data sets are not already ppp.objects")
		if(is.null(cases$x)||is.null(controls$y)||is.null(controls$x)||is.null(cases$y)) stop("data lists must both have two components named 'x' and 'y'")
		if((length(cases$x)!=length(cases$y))||(length(controls$x)!=length(controls$y))) stop("data components 'x' and 'y' are of unequal lengths in either 'cases' or 'controls'")
		if((length(cases$x)<10)||(length(controls$x)<10)) warning("less than 10 case or control observations!")
		cases.ppp <- ppp(x=cases$x,y=cases$y,window=WIN)
		controls.ppp <- ppp(x=controls$x,y=controls$y,window=WIN)
	} else if((class(cases)=="matrix")||(class(cases)=="data.frame")){
		if(is.null(WIN)) stop("must supply WIN if data sets are not already ppp.objects")
		if(length(na.omit(cases[,1]))!=length(na.omit(cases[,2]))) stop("data components 'x' and 'y' are of unequal lengths in 'cases'")
		if(length(na.omit(controls[,1]))!=length(na.omit(controls[,2]))) stop("data components 'x' and 'y' are of unequal lengths in 'controls'")
		if((nrow(cases)<10)||(nrow(controls)<10)) warning("less than 10 case or control observations!")
		cases.ppp <- ppp(x=cases[,1],y=cases[,2],window=WIN)
		controls.ppp <- ppp(x=controls[,1],y=controls[,2],window=WIN)
	} else {
		stop("invalid data object")
	}
	
	W <- cases$window
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

	if(is.null(hlim)){
		suppressWarnings(tppp <- ppp(x=c(cases$x,controls$x),y=c(cases$y,controls$y),window=W))
		md <- min(nndist(unique(tppp)))
		hlim <- c(md,max(md*50,min(diff(W$xrange),diff(W$yrange))/6))
	}
	
	if(comment) cat(paste("...searching for optimal common h within (",hlim[1],",",hlim[2],")...\n",sep=""))
	hopt <- suppressWarnings(optimise(LSCV.risk.single,interval=hlim,cases=cases.ppp,controls=controls.ppp,method=method[1],res=res,edge=edge,inside=selector)$minimum)
	return(hopt)
}
