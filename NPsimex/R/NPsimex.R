## Non-Parametric density estimation for data with measurement error using SIMEX
## Xiaofeng Wang
## wangx6@ccf.org
## Cleveland Clinic Foundation

simex.density <- function(W, msigma, x, from, to, n.user=128, n.lambda=50, lambda="SJ", span=8, adjust=1, na.rm = FALSE, ...){
	CheckValidity <- function(W,na.rm){
	    if (!is.numeric(W)) 
	      stop("argument 'W' must be numeric")
	    name <- deparse(substitute(W))
	    W <- as.vector(W)
	    W.na <- is.na(W)
	    if (any(W.na)) {
	      if (na.rm)  W <- W[!W.na]
	      else stop("'W' contains missing values")
	    }
	    W.finite <- is.finite(W)
	    if (any(!W.finite)) {
	      W <- W[W.finite]
	    }
	    nW <- length(W)
	    list(W=W,nW=nW,name=name, has.na=any(W.na));
	}
	Wout = CheckValidity(W);
	W=Wout$W; nW=Wout$nW; name=Wout$name; has.na=Wout$has.na
	if(length(msigma)!=1 | any(msigma<=0))
	  stop("'msigma' should a positive number in the case of homoscedastic error.");
	if (nW < 10) 
	  stop("The sample size is too small!")

	if(missing(x)){
	 if (missing(from)) 
	   from <- min(W)
	 if (missing(to)) 
	   to <- max(W)
	 if (!is.finite(from)) 
	   stop("non-finite 'from'")
	 if (!is.finite(to)) 
	   stop("non-finite 'to'")
	 if(from >= to){
	   stop("'from' is not smaller than 'to'!");
	 }else{
	   tt=seq(from, to, length=n.user);
	 }
	}else{
	   tt=sort(as.vector(x));
	}
	
	if (is.character(lambda)) {
        start <- switch(tolower(lambda), 
				nrd0 = lambda.select(W,msigma, method="nrd0"), 
				nrd = lambda.select(W,msigma, method="nrd"), 
            	ucv = lambda.select(W,msigma, method="ucv"), 
				bcv = lambda.select(W,msigma, method="bcv"), 
				sj = lambda.select(W,msigma, method="sj"), 
				`sj-ste` = lambda.select(W,msigma, method="sj-ste"), 
				`sj-dpi` = lambda.select(W,msigma, method="sj-dpi"), 
            	stop("unknown bandwidth rule"))
    } else {
		start <- lambda
		}
    if (!is.finite(start) | start <= 0) 
        stop("non-positive or non-finite 'lambda'!")
    
	start <- adjust * start
	end <- start+span
	lambda <- seq(start,end,length=n.lambda)
	# What <- matrix(0,nrow=length(tt),ncol=n.lambda)
	# ker <- function(x,t,b){mean(dnorm((t-x)/b)/b)}
	# for (i in 1:n.lambda) {
	# 	What[,i] <- apply(matrix(tt, ncol=1),1,ker,x=W, b=msigma*sqrt(lambda[i]))
	#  	}
	ksmooth <- function(x, msigma, tt, lambda){
		ll <- sqrt(lambda)
		n <- length(x)
		ntt <- length(tt)
		nll <- length(ll)
		nest <- ntt*nll
		est <- .C("SimexKern", as.double(x), as.integer(n), as.double(tt), as.integer(ntt), as.double(ll),  as.integer(nll), as.double(msigma), result = double(nest))
		est[["result"]]
	} 
	What <- t(matrix(ksmooth(x=W,msigma=msigma,tt=tt,lambda=lambda), nrow=n.user))
	
	simex1.est <- function(y,lambda){
		a <- lm(y~lambda+I(lambda^2))
		est <- a$coeff[1]-a$coeff[2]+a$coeff[3]
		return(est)
		}
	simex.est <- apply(as.matrix(What),2,simex1.est,lambda=lambda)
	simex.est[simex.est<0]  <- 0
	trapez.int <- function(x,y){
	  n <- length(x)
	  int <- (x[n]-x[1])/(n-1)*(sum(y)-(y[1]+y[n])/2)
	  return(int)
	}
	simex.int <- trapez.int(tt,simex.est)
	fhat <- simex.est/simex.int
	
	return(structure(list(x = tt, y = fhat, lambda1 = start, span=span, n.lambda=n.lambda, n = nW,
                        call = match.call(), data.name = name,
                        has.na = has.na), class = "NPsimex"))
}

simex.H.density <- function(W, msigma, x, from, to, n.user=128, n.lambda=50, lambda="SJ", span=8, adjust=1, na.rm = FALSE, ...){
	CheckValidity <- function(W,msigma,na.rm){
	    if (!is.numeric(W) | !is.numeric(msigma)) 
	      stop("argument 'W' or 'msigma' must be numeric.")
		if (length(W) != length(msigma)) stop("The lengths of 'W' and 'msigma' do not equal.")
	    name <- deparse(substitute(W))
	    W.msig <- cbind(W,msigma)
	    W.na <- is.na(W)
	    if (any(W.na)) {
	      if (na.rm)  W.msig <- W.msig[!W.na,]
	      else stop("'W' contains missing values")
	    }
	    msig.finite <- is.finite(W.msig[,2])
	    if (any(!msig.finite)) {
	      W.msig <- W.msig[msig.finite,]
	    }
	    nW <- nrow(W.msig)
	    list(W=W.msig[,1], msigma=W.msig[,2], nW=nW, name=name, has.na=any(W.na));
	}
	Wout = CheckValidity(W,msigma,na.rm);
	W=Wout$W; msigma=Wout$msigma; nW=Wout$nW; name=Wout$name; has.na=Wout$has.na
	
	if(any(msigma<=0)) stop("Standard deviations should be positive.")
	if (nW < 10) stop("The sample size is too small!")

	if(missing(x)){
	 if (missing(from)) 
	   from <- min(W)
	 if (missing(to)) 
	   to <- max(W)
	 if (!is.finite(from)) 
	   stop("non-finite 'from'")
	 if (!is.finite(to)) 
	   stop("non-finite 'to'")
	 if(from >= to){
	   stop("'from' is not smaller than 'to'!");
	 }else{
	   tt=seq(from, to, length=n.user);
	 }
	}else{
	   tt=sort(as.vector(x));
	}
	
	if (is.character(lambda)) {
        start <- switch(tolower(lambda), 
				nrd0 = lambda.select(W,msigma, method="nrd0"), 
				nrd = lambda.select(W,msigma, method="nrd"), 
            	ucv = lambda.select(W,msigma, method="ucv"), 
				bcv = lambda.select(W,msigma, method="bcv"), 
				sj = lambda.select(W,msigma, method="sj"), 
				`sj-ste` = lambda.select(W,msigma, method="sj-ste"), 
				`sj-dpi` = lambda.select(W,msigma, method="sj-dpi"), 
            	stop("unknown bandwidth rule"))
    } else {
		start <- lambda
		}
    if (!is.finite(start) | start <= 0) 
        stop("non-positive or non-finite 'lambda'!")
    
	start <- adjust * start
	end <- start+span
	lambda <- seq(start,end,length=n.lambda)
	# What <- matrix(0,nrow=length(tt),ncol=n.lambda)
	# ker <- function(x,t,b){mean(dnorm((t-x)/b)/b)}
	# for (i in 1:n.lambda) {
	# 	What[,i] <- apply(matrix(tt, ncol=1),1,ker,x=W, b=msigma*sqrt(lambda[i]))
	#  	}
	Hksmooth <- function(x, msigma, tt, lambda){
		ll <- sqrt(lambda)
		n <- length(x)
		ntt <- length(tt)
		nll <- length(ll)
		nest <- ntt*nll
		est <- .C("SimexKernH", as.double(x), as.integer(n), as.double(tt), as.integer(ntt), as.double(ll),  as.integer(nll), as.double(msigma), result = double(nest))
		est[["result"]]
	} 
	What <- t(matrix(Hksmooth(x=W,msigma=msigma,tt=tt,lambda=lambda), nrow=n.user))
	
	simex1.est <- function(y,lambda){
		a <- lm(y~lambda+I(lambda^2))
		est <- a$coeff[1]-a$coeff[2]+a$coeff[3]
		return(est)
		}
	simex.est <- apply(as.matrix(What),2,simex1.est,lambda=lambda)
	simex.est[simex.est<0]  <- 0
	trapez.int <- function(x,y){
	  n <- length(x)
	  int <- (x[n]-x[1])/(n-1)*(sum(y)-(y[1]+y[n])/2)
	  return(int)
	}
	simex.int <- trapez.int(tt,simex.est)
	fhat <- simex.est/simex.int
	
	return(structure(list(x = tt, y = fhat, lambda1 = start, span=span, n.lambda=n.lambda, n = nW,
                        call = match.call(), data.name = name,
                        has.na = has.na), class = "NPsimex"))
}

lambda.select <- function(W, msigma, method="SJ", na.rm = FALSE, ...){
    if (length(list(...)) > 0) 
        warning("non-matched further arguments are disregarded")
    if (!is.numeric(W)) 
        stop("argument 'W' must be numeric")
    name <- deparse(substitute(W))
    W <- as.vector(W)
    W.na <- is.na(W)
    if (any(W.na)) {
        if (na.rm) 
            W <- W[!W.na]
        else stop("'W' contains missing values")
    }
    N <- nW <- length(W)
    W.finite <- is.finite(W)
    if (any(!W.finite)) {
        W <- W[W.finite]
        nW <- length(W)
    }
    msigma.bar <- sqrt(mean(msigma^2))
    h0 <- switch(tolower(method), nrd0 = bw.nrd0(W), nrd = bw.nrd(W), 
        ucv = bw.ucv(W), bcv = bw.bcv(W), sj = , `sj-ste` = bw.SJ(W, 
            method = "ste"), `sj-dpi` = bw.SJ(W, method = "dpi"), 
        stop("unknown bandwidth rule!"))
	return((h0/msigma.bar)^2)
	}

print.NPsimex <- function (x, digits = NULL, ...){
    cat("\nCall:\n\t", deparse(x$call), "\n\nData: ", x$data.name, 
        " (", x$n, " obs.);", "\t lambda1 = ", formatC(x$lambda1, 
            digits = digits),"\t;", "\t span = ", formatC(x$span, 
	            digits = digits), "\n\n", sep = "")
    print(summary(as.data.frame(x[c("x", "y")])), digits = digits, 
        ...)
    invisible(x)
	}

plot.NPsimex <- function(x, type="l", xlab="x", ylab="y",  lwd=2, ...){
	plot.default(x, type=type, xlab=xlab, ylab=ylab, lwd=lwd,...)
	}

span.select <- function(W, msigma, span=c(2,4,6,8,10,12,25), from=min(W), to=max(W), n.user=128, n.lambda=50, lambda="SJ", bw="SJ", adjust=1, na.rm = FALSE, ...){
	if (length(msigma) != 1 || msigma <= 0) stop("'msigma' should be a positive numeric number.")
	nlength <- length(span)
	W.den <- density(W, bw=bw, from=from, to=to, n=n.user, ...)
	ISE <- rep(0, nlength)
	conden <- function(tt, ftt, msigma){
		conden <- apply(as.matrix(tt), 1, function(x) trapez.int(tt, ftt*dnorm(x-tt, sd=msigma)))
		return(conden)
	}
	trapez.int <- function(x,y){
		n <- length(x)
		int <- (x[n]-x[1])/(n-1)*(sum(y)-(y[1]+y[n])/2)
		return(int)
	}
	for (i in 1:nlength){
		X.simex <- simex.density(W, msigma=msigma, adjust=adjust, n.lambda=n.lambda, span=span[i], from=from, to=to, n.user=n.user, ...)
		W.con <- conden(X.simex$x, X.simex$y, msigma=msigma)
		ISE[i] <- trapez.int(X.simex$x, (W.con-W.den$y)^2)
	}
	return(list(span=span, ISE=ISE))
}


span.H.select <- function(W, msigma, span=c(2,4,6,8,10,12,25), approx=FALSE, from=min(W), to=max(W), n.user=128, n.lambda=50, lambda="SJ", bw="SJ", adjust=1, na.rm = FALSE, ...){
	if (length(W) != length(msigma)) stop("The lengths of 'W' and 'msigma' do not equal.")
	nlength <- length(span)
	W.den <- density(W, bw=bw, from=from, to=to, n=n.user, ...)
	ISE <- rep(0, nlength)
	conden <- function(tt, ftt, msigma){
		conden <- apply(as.matrix(tt), 1, function(x) trapez.int(tt, ftt*dnorm(x-tt, sd=msigma)))
		return(conden)
	}
	trapez.int <- function(x,y){
		n <- length(x)
		int <- (x[n]-x[1])/(n-1)*(sum(y)-(y[1]+y[n])/2)
		return(int)
	}
	if (approx) {
		msigma.bar <- sqrt(mean(msigma^2))
		for (i in 1:nlength){
			X.simex <- simex.H.density(W, msigma=msigma, adjust=adjust, n.lambda=n.lambda, span=span[i], from=from, to=to, n.user=n.user, ...)
			W.con <- conden(X.simex$x, X.simex$y, msigma=msigma.bar)
			ISE[i] <- trapez.int(X.simex$x, (W.con-W.den$y)^2)
		}
	} else {
		msigma1 <- unique(msigma)
		for (i in 1:nlength){
			X.simex <- simex.H.density(W, msigma=msigma, adjust=adjust, n.lambda=n.lambda, span=span[i], from=from, to=to, n.user=n.user, ...)
			W.con <- rowMeans(apply(as.matrix(msigma1),1, conden, tt=X.simex$x, ftt=X.simex$y))
			ISE[i] <- trapez.int(X.simex$x, (W.con-W.den$y)^2)
		}
	}
	return(list(span=span, ISE=ISE))
}


