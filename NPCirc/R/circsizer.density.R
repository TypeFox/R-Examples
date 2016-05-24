circsizer.density <- function(x,bws,ngrid=250,alpha=0.05,B=500,log.scale=TRUE,display=TRUE){
	data <- x
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
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
	x.na <- is.na(x)
	if (sum(x.na)>0) warning("Missing values were removed")
	x <- x[!x.na]
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (!is.numeric(bws)) stop("Argument 'bws' must be numeric")
	if (length(bws)==1) stop("Length of argument 'bws' must be at least 2")
	if (any(bws<0)) warning("Values of vector 'bws' must be positive. Values smaller than 0 were removed")
	bws <- sort(bws[bws>=0])
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

	
	t <- seq(0,2*pi,length=ngrid)
	t_x <- outer(t,x,"-")
	nbws <- length(bws.plot) 
	COLOR<-matrix("magenta",nbws,ngrid)
	LL <- UL <- ESS <- matrix(NA,nbws,ngrid)
	start <- 1
	if ((!log.scale & bws[1]==0) | (log.scale & bws[nbws]==0)) {
		LL[1,] <- UL[1,] <-0
		ESS[1] <- n
		start <- 2
	}
	for (i in start:nbws){

		bw <- bws[i]
		exptx <- exp(cos(t_x)*bw)
		exptxsin <- exptx*sin(t_x)  
		pibessel <- 2*pi*besselI(bw,0)

		# Derivative estimate
		deriv_fhat <- -bw*rowMeans(exptxsin)/pibessel 

		# Derivative of the von Mises kernel
		deriv_kernel <- -bw*exptxsin/pibessel

		if (sum(is.na(deriv_fhat))>0 | sum(is.na(deriv_kernel))>0 | sum(deriv_kernel==Inf)>0){
			stop("Values of the smoothing parameter too large")
		}

		# Standard deviation	
		sd <- sqrt((rowMeans(deriv_kernel^2)-rowMeans(deriv_kernel)^2)/n)

		# Effective Sample Size
		ess <-rowSums(exptx)/exp(bw) 
		
		# Bootstrap-t confidence interval	
		zstar<-matrix(NA,B,ngrid)
		for (b in 1:B){
			xstar<- sample(x,replace=T)
			t_xstar <- outer(t,xstar,"-")
			exptxstar <- exp(cos(t_xstar)*bw)
			exptxstarsin <- exptxstar*sin(t_xstar)
			deriv_fhat_xstar <- -bw*rowMeans(exptxstarsin)/pibessel
			deriv_kernel <- -bw*exptxstarsin/pibessel 
			sdest <- sqrt((rowMeans(deriv_kernel^2)-rowMeans(deriv_kernel)^2)/n)
			zstar[b,]<-(deriv_fhat_xstar-deriv_fhat)/sdest
		}
		quant <- apply(zstar, 2, function(z) quantile(z, probs=c(alpha/2,1-alpha/2),type=1))
		ll <- deriv_fhat-quant[2,]*sd
		ul <- deriv_fhat-quant[1,]*sd
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
