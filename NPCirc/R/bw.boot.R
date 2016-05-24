bw.boot <- function(x,lower=0,upper=100,np=500,tol=0.1){
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	attr(x, "class") <- attr(x, "circularp") <- NULL
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	x.na <- is.na(x)
	if (sum(x.na)>0) warning("Missing values were removed")
	x <- x[!x.na]
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (!is.numeric(upper)){ 
		warning("argument 'upper' must be numeric. Default upper boundary was used")
		upper <- 100
	}
	if (!is.numeric(lower)){
		warning("argument 'lower' must be numeric. Default lower boundary was used")
		lower <- 0
	}
	if (lower<0 | lower>=upper){
      	warning("The boundaries must be positive and 'lower' must be smaller that 'upper'. Default boundaries were used")
		upper <- 100
		lower <- 0
	}
	if (!is.numeric(np)) stop("argument 'np' must be numeric")
	if (np<=0){ 
		warning("'np' must be positive. Default value of 'np' was used")
		np <- 500
	}
	t <- seq(0,2*pi,length=np)
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")
	MISEboot<-function(t,x,bw){
		n <- length(x)
		costx <- cos(outer(t,x,"-"))
		exptx <- exp(costx*bw)
		fhat <- apply(exptx,1,sum)/(n*2*pi*besselI(bw,0))
		EBf <- apply((besselI(abs(2*bw*cos(outer(t,x,"-")/2)),0)/besselI(bw,0)),1,sum)/(2*pi*n*besselI(bw,0))
		EBff2 <-(apply(besselI(bw*sqrt(5+4*costx),0),1,sum)/besselI(bw,0))/((2*pi*n*besselI(bw,0))^2)+(EBf-fhat)^2 - EBf^2/n
		return(int.Simp(t,EBff2))
	}
	bw<-optimize(function(bw) MISEboot(t,x,bw),interval=c(lower,upper),tol=tol)$minimum
	return(bw)
}
