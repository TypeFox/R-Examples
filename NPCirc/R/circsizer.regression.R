circsizer.regression <- function(x,y,bws=NULL,adjust=2,ngrid=150,alpha=0.05,B=500,B2=250,log.scale=TRUE,display=TRUE){
	data <- cbind(x,y)
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.numeric(y)) stop("argument 'y' must be numeric")
	if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
	xcircularp <- attr(x, "circularp")
	if (is.null(xcircularp)){
		type <- 3
		zero <- 0
		clockwise <- FALSE
	}else{
        	units <- xcircularp$units
		zero <- xcircularp$zero
		rotation <- xcircularp$rotation
		clockwise <- ifelse(rotation=="counter", FALSE, TRUE)
        	template <- xcircularp$template
    		if (template == "geographics") type <- 1
		else if (template == "clock24") type <- 2
    		else if (units == "radians") type <- 3
		else if (units == "degrees") type <- 4
    	}
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	attr(x, "class") <- attr(x, "circularp") <- NULL
	nax <- is.na(x)
	nay <- is.na(y)
	x <- x[!nax & !nay]
	y <- y[!nax & !nay]
	if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.")
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (is.null(bws)){
		bw <- bw.reg.circ.lin(x,y,method="LL")
		bws<-c(bw/adjust,bw*adjust)
	}else{
		if (!is.numeric(bws)) stop("argument 'bws' must be numeric")
		if (length(bws)==1) stop("Length of argument 'bws' must be at least 2")
		if (any(bws<0)) warning("Values of vector 'bws' must be positive. Values smaller than 0 were removed")
		bws <- sort(bws[bws>=0]) 
	}
	if (!log.scale){ 
		bws.eq <- seq(min(bws),max(bws),length=length(bws))
		if (all(bws==bws.eq)==FALSE){
			bws <- bws.eq
			warning ("Values of 'bws' are not equally spaced. CircSiZer is computed for a equally spaced grid between 
			the minimum and maximum values of 'bws'")
		}
		bws.plot <- bws
	}else{ 
		bws.plot <- seq(-log10(max(bws[bws!=0])),-log10(min(bws[bws!=0])),length=length(bws))
		bws <- 10^(-bws.plot)
		if(length(unique(bws))==1){
			warning("Argument 'bws' is not valid. A default grid was used")
			bws.plot<-seq(-2,1,by=0.1)
			bws<-10^(-bws.plot)
		}
	}
	if (!is.numeric(ngrid)) stop("argument 'ngrid' must be numeric")
	if (ngrid<=0){
		warning("'ngrid' must be positive integer. Default value of 'ngrid' was used")
		ngrid=250
	}
	if (!is.numeric(alpha)) stop("argument 'alpha' must be numeric")
	if (alpha<0 | alpha>1){
		warning("'alpha' must be in the interval [0,1]. Default value of 'alpha' was used")
		alpha=0.05
	}
	if (!is.numeric(B)) stop("argument 'B' must be numeric")
	if (B<=0){
		warning("'B' must be positive. Default value of 'B' was used")
		B=500
	}
	if (!is.numeric(B2)) stop("argument 'B' must be numeric")
	if (B2<=0){
		warning("'B2' must be positive. Default value of 'B2' was used")
		B2=250
	}

	# Auxiliary function that computes the estimator of the derivate of the regression function
	der.cl<-function(x,y,bw){
		x_t <- matrix(x,nrow=n)%*%v1t - v1xt
		m<-sin(x_t)
		kxt<-exp(bw*cos(x_t))
		sn1<-colSums(kxt*m)   
		sn2<-colSums(kxt*m^2) 
		sn3<-colSums(kxt)     
		titawy1<-colSums(kxt*y)   
		titawy2<-colSums(kxt*m*y)
		beta1<- -sn1*titawy1+sn3*titawy2
		det<-sn3*sn2-sn1^2 
		return(beta1/det)
	}

	t <- seq(0,2*pi,length=ngrid)
	v1t<-rep(1,ngrid)
	v1xt<-matrix(rep(1,n),nrow=n)%*%t
	nbws <- length(bws.plot) 
	COLOR <- matrix("magenta",nbws,ngrid)
	LL <- UL <- ESS <- matrix(NA,nbws,ngrid)
	start <- 1
	if ((!log.scale & bws[1]==0) | (log.scale & bws[nbws]==0)){
		LL[1,]<-UL[1,]<-0
		ESS[1] <- n
		start <- 2
	}
	for (i in start:nbws){
		bw <- bws[i]
		# Estimate of the derivate of the regression function
		fhat <- der.cl(x,y,bw) 
		if (sum(is.na(fhat))>0){
			stop("Values of the smoothing parameter too large")
		}
		# Standard deviation
		Xstar<-Ystar<-matrix(0,B,n)
		der<-matrix(0,B,ngrid)
		ind<-matrix(sample(n,n*B,replace=T),n,B)
		for (k in 1:B){
			indk<-ind[,k]
			xstar<-x[indk]
			ystar<-y[indk]
			der[k,]<-der.cl(xstar,ystar,bw) 
			Xstar[k,]<-xstar
			Ystar[k,]<-ystar
		}
		sd<-sqrt(colMeans(der^2)-colMeans(der)^2)
		# Effective Sample Size
		exptx<-exp(bw*cos(outer(t,x,"-")))
		ess<-rowSums(exptx)/exp(bw)
		# Bootstrap-t confidence interval		
		zstar<-matrix(0,B,ngrid)
		derstarstar<-matrix(0,B2,ngrid)
		for (l in 1:B){
			ind<-matrix(sample(n,n*B2,replace=T),n,B2)
			for (k in 1:B2){
				indk<-ind[,k]
				xstarstar<-Xstar[l,][indk]
				ystarstar<-Ystar[l,][indk]
				derstarstar[k,] <- der.cl(xstarstar,ystarstar,bw)
			}
			sdest<-sqrt(colMeans(derstarstar^2)-colMeans(derstarstar)^2)
			zstar[l,]<-(der[l,]-fhat)/sdest
		}
		quant <- apply(zstar, 2, function(z) quantile(z, probs=c(alpha/2,1-alpha/2),type=1))
		ll <- fhat-quant[2,]*sd
		ul <- fhat-quant[1,]*sd
		LL[i,] <- ll
		UL[i,] <- ul
		ESS[i,] <- ess
		COLOR[i,ll>0] <- "darkblue"
		COLOR[i,ul<0] <- "red"
		COLOR[i,ess<5] <- "grey"
	}
	circsizer.object <- structure(list(data=data, ngrid=ngrid, bw=bws.plot, log.scale=log.scale, 
	CI=list(LL,UL,ESS), col=COLOR, call = match.call(), data.name = "circsizer.object"), class = "circsizer")
	if (display) circsizer.map(circsizer.object, type=type, zero=zero, clockwise=clockwise)
	return(circsizer.object)
}