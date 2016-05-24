bw.CV<-function(x,method="LCV",lower=0,upper=50,tol=1e-2,np=500){
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
		upper <- 50
	}
	if (!is.numeric(lower)){
		warning("argument 'lower' must be numeric. Default lower boundary was used")
		lower <- 0
	}
	if (lower<0 | lower>=upper){
      	warning("The boundaries must be positive and 'lower' must be smaller that 'upper'. Default boundaries were used")
		upper <- 50
		lower <- 0
	}
	if (!is.numeric(np)) stop("argument 'np' must be numeric")
	if (!is.numeric(tol)) stop("argument 'tol' must be numeric")
	if (method=="LSCV"){
		if (np<=0){ 
			warning("'np' must be positive. Default value of 'np' was used")
			np <- 500
		}
		t<-seq(0,2*pi,length=np)
		lscv<-function(x,bw){
			est<-DensCircRad(x,t,bw)
			n<-length(x)
			cv<-numeric(n)
			for (j in 1:n){
				cv[j]<-DensCircRad(x[-j],x[j],bw)
			}
			return(int.Simp(t,est^2)-2*mean(cv))	
		}
		bw <- optimize(function(h)lscv(x,h),interval=c(lower,upper),tol=tol)$minimum
	}else{
		cv <- function(x,bw){ 
			n<-length(x)
			logdens<-numeric(n)
			for (j in 1:n){
				logdens[j]<-log(DensCircRad(x[-j],x[j],bw))
			}
			return(sum(logdens))
		}
		bw <- optimize(function(h)cv(x,h),interval=c(lower,upper),tol=tol,maximum=TRUE)$maximum
	}
	if (bw < lower + tol | bw > upper - tol) 
      warning("minimum/maximum occurred at one end of the range")
	return(bw)
}

