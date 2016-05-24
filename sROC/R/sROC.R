## Non-Parametric Smooth ROC Curves for Continuous Data
## Xiao-Feng Wang
## wangx6@ccf.org
## Cleveland Clinic Lerner Research Institute

## kernel estimator for cumulative distribution function (CDF)

kCDF <- function(x, bw="pi_ucv", adjust=1, kernel=c("normal", "epanechnikov"), xgrid, ngrid=256, from, to, cut=3, na.rm = FALSE, ...){
	if (length(list(...))) 
        warning("non-matched further arguments are disregarded")
	kernel <- match.arg(kernel)
    CheckValidity <- function(x,na.rm){
	    if (!is.numeric(x)) 
	      stop("argument 'x' must be numeric")
	    name <- deparse(substitute(x))
	    x <- as.vector(x)
	    x.na <- is.na(x)
	    if (any(x.na)) {
	      if (na.rm)  x <- x[!x.na]
	      else stop("'x' contains missing values")
	    }
	    x.finite <- is.finite(x)
	    if (any(!x.finite)) {
	      x <- x[x.finite]
	    }
	    n<- length(x)
	    list(x=x,n=n, name=name, has.na=any(x.na));
	}
	data <- x
	xout = CheckValidity(x, na.rm=na.rm)
	x=xout$x; n=xout$n; name=xout$name; has.na=xout$has.na
	if (n < 5) stop("The sample size is too small!")
	
	if (is.character(bw)) {
        bw <- switch(tolower(bw), 
			ncdf = bw.CDF(x, "ncdf"), 
			npdf = bw.CDF(x, "npdf"), 
			pi_nrd0 = bw.CDF.pi(x, pilot="nrd0"),
			pi_nrd = bw.CDF.pi(x, pilot="nrd"),
			pi_ucv = bw.CDF.pi(x, pilot="ucv"),
			pi_bcv = bw.CDF.pi(x, pilot="bcv"),
			pi_sj = bw.CDF.pi(x, pilot="sj"),
			pi_onestage = bw.CDF.pi(x, pilot="onestage"),
            stop("unknown bandwidth rule"))
    }
    if (!is.finite(bw)) 
        stop("non-finite 'bw'")
    bw <- adjust * bw
	
	if(missing(xgrid)){
	 if (missing(from)) 
	   from <- min(x) - cut*bw
	 if (missing(to)) 
	   to <- max(x) + cut*bw
	 if (!is.finite(from)) 
	   stop("non-finite 'from'")
	 if (!is.finite(to)) 
	   stop("non-finite 'to'")
	 if(from >= to){
	   stop("'from' is not smaller than 'to'!")
	 }else{
	   xgrid=seq(from, to, length=ngrid);
	 }
	}else{
	   xgrid=sort(as.vector(xgrid));
	}
	ngrid <- length(xgrid)

	kCDF <- switch(substr(tolower(kernel), 1, 4), 
		norm = .C("NKern", x=as.double(x), n=as.integer(n), xgrid=as.double(xgrid), ngrid=as.integer(ngrid), bw=as.double(bw), Fhat=double(ngrid)),
		epan = .C("EKern", x=as.double(x), n=as.integer(n), xgrid=as.double(xgrid), ngrid=as.integer(ngrid), bw=as.double(bw), Fhat=double(ngrid)),
		stop("The specified kernel type is not supported!"))

	return(structure(list(x = kCDF[["xgrid"]], Fhat = kCDF[["Fhat"]], bw = bw, 
        n = n, call = match.call(), data.name = name, data=data, has.na=has.na), 
        class = "CDF"))
}

print.CDF <- function (x, digits = NULL, ...) 
{
    cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
        " (", x$n, " obs.);", "\tBandwidth 'bw' = ", formatC(x$bw, 
            digits = digits), "\n\n", sep = "")
    print(summary(as.data.frame(x[c("x", "Fhat")])), digits = digits, ...)
    invisible(x)
}

## Bandwidth selection for kCDF

bw.CDF <- function (x, method="npdf") 
{
	n <- length(x)
    if (n < 2L) 
        stop("need at least 2 data points")
	xsd <- sd(x)
	bw <- switch(tolower(method), ncdf=(180*sqrt(pi)/(7*n))^(1/3)*min(xsd, IQR(x)/1.349),
			npdf=(4/(3*n))^(1/5)*min(xsd, IQR(x)/1.349),
			stop("Unknown bandwidth rule!!!"))
    return(bw)
}

bw.CDF.pi <- function (x, pilot="UCV") 
{
	if (length(x) < 2L) stop("need at least 2 data points")
	n <- length(x)
	xsd <- sd(x)
    xgrid <- x
    ngrid <- length(xgrid)
	g <- switch(tolower(pilot), 
			nrd0=bw.nrd0(x),
			nrd=bw.nrd(x),
			ucv=bw.ucv(x),
			bcv=bw.bcv(x),
			sj=bw.SJ(x),
			onestage=(3/(n*3/(8*sqrt(pi))*(min(xsd, IQR(x)/1.349))^(-5)))^(1/5), 
			stop("Unknown bandwidth rule!!!"))
	phi <- .C("phi2", x=as.double(x), n=as.integer(n), xgrid=as.double(xgrid), ngrid=as.integer(ngrid), bw=as.double(g), phi=double(ngrid))
    phi2 <- mean(phi[["phi"]])
	bw <- (45/(-7*n*phi2))^(1/3)
    return(bw)
}

## ----------------Pointwise Confidence Intervals for Smooth CDF --------------
## confidence level 100(1 − \alpha) %

CI.CDF <- function(CDF, alpha=0.05){
	x <- CDF$x
	Fhat <- CDF$Fhat
	n <- CDF$n
	Fhat.SD <- sqrt(Fhat*(1-Fhat)/n)
	crit <- qnorm(1-alpha/2)
	Fhat.upper <- Fhat + crit*Fhat.SD
	Fhat.lower <- Fhat - crit*Fhat.SD
	return(list(x=x, Fhat=Fhat, Fhat.upper=Fhat.upper, Fhat.lower=Fhat.lower, alpha=alpha))
}

plot.CDF  <- function (x, CI=TRUE, alpha=0.05, main = NULL, xlab = NULL, ylab = "CDF", lwd=2, lty=1, ...) 
{
    if (is.null(xlab)) 
        xlab <- paste("N =", x$n, "  Bandwidth =", formatC(x$bw))
    if (is.null(main)) 
        main <- deparse(x$call)
	if (CI) {
		result <- CI.CDF(x, alpha=alpha)
		plot.default(result$x, result$Fhat, main = main, xlab = xlab, ylim= c(0,1), ylab = ylab, type = "n", lwd=lwd, ...)
		lines(result$x, result$Fhat.lower, lwd=lwd, lty=lty, ...)
		lines(result$x, result$Fhat.upper, lwd=lwd, lty=lty, ...)
		polygon(c(result$x, rev(result$x)), c(result$Fhat.lower, 
		     rev(result$Fhat.upper)), col = "grey", 
		     border = FALSE)
		lines(result$x, result$Fhat, lwd=lwd, lty=lty, ...)
	} else {
		plot.default(x$x, x$Fhat, main = main, xlab = xlab, ylim= c(0,1), ylab = ylab, type = "l", lwd=lwd, ...)
	}
    invisible(NULL)
}



## kernel estimator for receiver operating characteristic (ROC) Curve

kROC <- function(x, y, bw.x="pi_ucv", bw.y="pi_ucv", adjust=1, kernel=c("normal", "epanechnikov"), xgrid, ngrid=256, from, to, cut=3, na.rm = FALSE, ...){
	if (length(list(...))) 
        warning("non-matched further arguments are disregarded")

	kernel <- match.arg(kernel)
    CheckValidity <- function(x,na.rm){
	    if (!is.numeric(x)) 
	      stop("argument 'x' or 'y' must be numeric")
	    name <- deparse(substitute(x))
	    x <- as.vector(x)
	    x.na <- is.na(x)
	    if (any(x.na)) {
	      if (na.rm)  x <- x[!x.na]
	      else stop("'x' or 'y' contains missing values")
	    }
	    x.finite <- is.finite(x)
	    if (any(!x.finite)) {
	      x <- x[x.finite]
	    }
	    n <- length(x)
	    list(x=x,n=n, name=name, has.na=any(x.na));
	}
	xout = CheckValidity(x, na.rm=na.rm)
	x=xout$x; nx=xout$n; x.name=xout$name; x.has.na=xout$has.na
	yout = CheckValidity(y, na.rm=na.rm)
	y=yout$x; ny=yout$n; y.name=yout$name; y.has.na=yout$has.na
	if (nx < 5 | ny < 5) stop("The sample size of x or y is too small!")

	if (is.character(bw.x)) {
        bw.x <- switch(tolower(bw.x), 
			ncdf = bw.CDF(x, "ncdf"), 
			npdf = bw.CDF(x, "npdf"), 
			pi_nrd0 = bw.CDF.pi(x, pilot="nrd0"),
			pi_nrd = bw.CDF.pi(x, pilot="nrd"),
			pi_ucv = bw.CDF.pi(x, pilot="ucv"),
			pi_bcv = bw.CDF.pi(x, pilot="bcv"),
			pi_sj = bw.CDF.pi(x, pilot="sj"),
			pi_onestage = bw.CDF.pi(x, pilot="onestage"),
            stop("unknown bandwidth rule"))
    }
	if (is.character(bw.y)) {
        bw.y <- switch(tolower(bw.y), 
			ncdf = bw.CDF(y, "ncdf"), 
			npdf = bw.CDF(y, "npdf"), 
			pi_nrd0 = bw.CDF.pi(y, pilot="nrd0"),
			pi_nrd = bw.CDF.pi(y, pilot="nrd"),
			pi_ucv = bw.CDF.pi(y, pilot="ucv"),
			pi_bcv = bw.CDF.pi(y, pilot="bcv"),
			pi_sj = bw.CDF.pi(y, pilot="sj"),
			pi_onestage = bw.CDF.pi(y, pilot="onestage"),
            stop("unknown bandwidth rule"))
    }

	if(missing(xgrid)){
	 if (missing(from)) 
	   from <- min(c(x,y)) - cut*min(c(bw.x,bw.y))
	 if (missing(to)) 
	   to <- max(c(x,y)) + cut*min(c(bw.x,bw.y))
	 if (!is.finite(from)) 
	   stop("non-finite 'from'")
	 if (!is.finite(to)) 
	   stop("non-finite 'to'")
	 if(from >= to){
	   stop("'from' is not smaller than 'to'!")
	 }else{
	   xgrid=seq(from, to, length=ngrid);
	 }
	}else{
	   xgrid=sort(as.vector(xgrid));
	}

	x.CDF <- kCDF(x,bw=bw.x, adjust=adjust, kernel=kernel, xgrid=xgrid)
	y.CDF <- kCDF(y,bw=bw.y, adjust=adjust, kernel=kernel, xgrid=xgrid)
	return(structure(list(FPR =1 - y.CDF$Fhat, TPR=1 - x.CDF$Fhat, bw.x = bw.x, bw.y = bw.y, 
        nx = nx, ny = ny, call = match.call(), x.data.name = x.name, y.data.name = y.name, x.has.na=x.has.na, y.has.na=y.has.na), 
        class = "ROC"))
}

print.ROC <- function (x, digits = NULL, ...) 
{
    cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$x.data.name, 
        " (", x$nx, " obs.) and ", x$y.data.name, " (", x$ny, " obs.);", "\n\n", sep = "")
	cat("\tBandwidth 'bw.x' = ", formatC(x$bw.x, digits = digits), " 'bw.y' = ", formatC(x$bw.y, 
		            digits = digits),"\n\n", sep = "")
    print(summary(as.data.frame(x[c("FPR", "TPR")])), digits = digits, 
        ...)
    invisible(x)
}

plot.ROC  <- function (x, main = NULL, diagonal = TRUE, xlab = "FPR", ylab = "TPR", type = "l", lwd=2, ...) 
{
    if (is.null(main)) 
        main <- "Smooth ROC curve"
    plot.default(x$FPR, x$TPR, main = main, xlab = xlab, ylab = ylab, xlim=c(0,1), ylim= c(0,1), type = type, lwd=lwd,
        ...)
	if (diagonal) abline(a=0,b=1)
    invisible(NULL)
}

AUC <- function(ROC, method="Simpson", ngrid=256){
	FPR <- ROC$FPR
	TPR <- ROC$TPR
	
	TrapezInt <- function(x, y, ngrid = 501, equal=FALSE){
		n <- length(x)
		if (n != length(y)) stop("'x' and 'y' have different length.")
		if (equal) {
			int <- (x[2]-x[1])*(sum(y)-y[1]/2-y[n]/2)
		} else {
	    	ap <- approx(x, y, n=ngrid)
	 		int <- (ap$x[2]-ap$x[1])*(sum(ap$y)-ap$y[1]/2-ap$y[ngrid]/2)
		}		
		return(int)
	}

	SimpsonInt <- function (x, y, ngrid = 256) 
		{
	    	if (length(x) != length(y)) stop("'x' and 'y' have different length.")
	    	ap <- approx(x, y, n = 2 * ngrid + 1)
	    	int <- (ap$x[2]-ap$x[1]) * (ap$y[2 * (1:ngrid) - 1] + 4 * ap$y[2 * (1:ngrid)] + ap$y[2 * (1:ngrid) + 1])/3
	    	return(sum(int))
		}

	AUC <- switch(substr(tolower(method), 1, 4),
	 			trap = TrapezInt(FPR,TPR,ngrid=ngrid),
	 			simp = SimpsonInt(FPR,TPR,ngrid=ngrid),
	            stop("unknown numerical integration method!"))
	return(structure(list(AUC=AUC), class="AUC"))
}

print.AUC <- function(x, digits = NULL, ...){
	cat("\tAUC for the ROC curve is ", formatC(x$AUC, digits = digits),".\n\n", sep = "")
}

