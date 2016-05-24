kern.reg.lin.circ<-function(x,y,t=NULL,bw=NULL,method="LL",len=250){
	name <- deparse(substitute(x))
	datax <- x
	datay <- y
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.numeric(y)) stop("argument 'y' must be numeric")
	if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
	if (is.circular(y)){
		datacircularp <- circularp(y)	
    	} else {
		datacircularp <- list(type = "angles", units = "radians", template = "none", modulo = "2pi", zero = 0, rotation = "counter")
    	}
	dc <- list()
      dc$type <- datacircularp$type
      dc$units <- datacircularp$units
	dc$template <- datacircularp$template
	dc$modulo <- datacircularp$modulo
	dc$zero <- datacircularp$zero
	dc$rotation <- datacircularp$rotation
	if (is.null(t)){
		t <- seq(min(x), max(x), length=len)
	}else{
		if (!is.numeric(t)) stop("argument 't' must be numeric")
		t.na <- is.na(t)
		t <- t[!t.na]
		if (sum(t.na)>0) warning("'t' contains missing values. They were removed")
	}
	y <- conversion.circular(y, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	nax <- is.na(x)
	nay <- is.na(y)
	x<-x[!nax & !nay]
	y<-y[!nax & !nay]
	if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.", "\n")
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (is.null(bw)){
		bw <- bw.reg.lin.circ(x,y,method)		
	}else{
		if (is.numeric(bw)){
			if (bw<0){
				warning("Argument 'bw' must be positive. The value of 'bw' was computed by using the plug--in rule")
				bw <- bw.reg.lin.circ(x,y,method)	
			}
		}else{
 			warning("Argument 'bw' must be numeric. The value of 'bw' was computed by using the plug--in rule")
			bw <- bw.reg.lin.circ(x,y,method)	
		}
	}
	attr(y, "class") <- attr(y, "circularp") <- NULL
	if (method=="NW") fhat <- RegLinCirc(x,y,t,bw,method="NW")
	else fhat <- RegLinCirc(x,y,t,bw,method="LL")
	fhat <- conversion.circular(circular(fhat), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
	structure(list(datax = datax, datay = datay, x = t, y = fhat, bw = bw, n = n, kernel = "gaussian", 
	call = match.call(), data.name = name, has.na = FALSE), class = "regression.circular")
}



