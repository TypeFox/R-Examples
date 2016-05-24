kern.reg.circ.circ<-function(x,y,t=NULL,bw=NULL,method="LL",from=circular(0),to=circular(2*pi),len=250){
	name <- deparse(substitute(x))
	datax <- x
	datay <- y
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.numeric(y)) stop("argument 'y' must be numeric")
	if (length(x) != length(y)) stop("'x' and 'y' must have the same number of observations")
	if (!is.null(t) && is.circular(t)) {
		datacircularp <- circularp(t)
	} else if (is.circular(x)){
		datacircularp <- circularp(x)
	} else if (is.circular(y)){
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
		if (!is.numeric(from)) stop("argument 'from' must be numeric")
		if (!is.numeric(to)) stop("argument 'to' must be numeric")
		if (!is.finite(from)) stop("non-finite 'from'")
		if (!is.finite(to)) stop("non-finite 'to'")
		if (!is.numeric(len)) stop("argument 'len' must be numeric")
		if (len <= 0) stop("argument 'len' must be integer and positive")
		from <- conversion.circular(from, units = "radians", zero = 0, rotation = "counter")
		attr(from, "class") <- attr(from, "circularp") <- NULL
		to <- conversion.circular(to, units = "radians", zero = 0, rotation = "counter")
    		attr(to, "class") <- attr(to, "circularp") <- NULL
		if (from>to) stop("argument 'from' must be smaller than argument 'to'")
		t <- circular(seq(from=from,to=to,length=len))
	}else{
		if (!is.numeric(t)) stop("argument 't' must be numeric")
		t.na <- is.na(t)
		t <- t[!t.na]
		if (sum(t.na)>0) warning("'t' contains missing values. They were removed")
	}
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	y <- conversion.circular(y, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	nax <- is.na(x)
	nay <- is.na(y)
	x<-x[!nax & !nay]
	y<-y[!nax & !nay]
	if ((sum(nax)+sum(nay))>0) warning("Missing values were removed.", "\n")
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
	if (is.null(bw)){
		bw <- bw.reg.circ.circ(x,y,method)		
	}else{
		if (is.numeric(bw)){
			if (bw<0){
				warning("Argument 'bw' must be positive. The value of 'bw' was computed by using the plug--in rule")
				bw <- bw.reg.circ.circ(x,y,method)	
			}
		}else{
 			warning("Argument 'bw' must be numeric. The value of 'bw' was computed by using the plug--in rule")
			bw <- bw.reg.circ.circ(x,y,method)	
		}
	}
	attr(x, "class") <- attr(x, "circularp") <- NULL
	attr(y, "class") <- attr(y, "circularp") <- NULL
	tt <- conversion.circular(t, dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
	t <- conversion.circular(t, units = "radians", modulo="2pi", zero = 0, rotation = "counter")
	attr(t, "class") <- attr(t, "circularp") <- NULL
	if (method=="NW") fhat <- RegCircCirc(x, y, t, bw, method="NW")
	else fhat <- RegCircCirc(x, y, t, bw, method="LL")
	fhat <- conversion.circular(circular(fhat), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
	structure(list(datax = datax, datay = datay, x = tt, y = fhat, bw = bw, n = n, kernel = "vonmises", 
	call = match.call(), data.name = name, has.na = FALSE), class = "regression.circular")
}



