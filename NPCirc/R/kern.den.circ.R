kern.den.circ<-function(x,t=NULL,bw=NULL,from=circular(0),to=circular(2*pi),len=250){
	name <- deparse(substitute(x))
	data <- x
	if (!is.numeric(x)) stop("argument 'x' must be numeric")
	if (!is.null(t) && is.circular(t)) {
		datacircularp <- circularp(t)
	} else if (is.circular(x)) {
		datacircularp <- circularp(x)
    	} else {
		datacircularp <- list(type = "angles", units = "radians", template = "none", modulo = "2pi", zero = 0, rotation = "counter")
    	}
	x <- conversion.circular(x, units = "radians", zero = 0, rotation = "counter", modulo = "2pi")
	if (is.null(bw)){
		bw<-bw.pi(x)		
	}else{
		if (is.numeric(bw)){
			if (bw<0){
				warning("Argument 'bw' must be positive. The value of 'bw' was computed by using the plug--in rule")
				bw<-bw.pi(x)
			}
		}else{
 			warning("argument 'bw' must be numeric. The value of 'bw' was computed by using the plug--in rule")
			bw<-bw.pi(x)
		}
	}
	attr(x, "class") <- attr(x, "circularp") <- NULL
	x.na <- is.na(x)
	if (sum(x.na)>0) warning("Missing values were removed")
	x <- x[!x.na]
	n <- length(x)
	if (n==0) stop("No observations (at least after removing missing values)")
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
	tt <- conversion.circular(t, dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
	t <- conversion.circular(t, units = "radians", modulo="2pi", zero = 0, rotation = "counter")
	attr(t, "class") <- attr(t, "circularp") <- NULL
	y <- DensCircRad(x, t, bw)
	structure(list(data = data, x = tt, y = y, bw = bw, n = n, kernel = "vonmises", 
	call = match.call(), data.name = name, has.na = FALSE), class = "density.circular")
}

