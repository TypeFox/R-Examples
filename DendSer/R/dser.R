



dser <- function(x,ser_weight,cost=costBAR, ...)
         UseMethod("dser")
         


dser.data.frame <- function(x,ser_weight,cost=costBAR,...) {
	dser.matrix(as.matrix(x),ser_weight,cost=cost,...)
	}
	
	
dser.matrix <-function(x,ser_weight,cost=costBAR,scale=TRUE,dmethod="euclidean",...) {
	if (!is.matrix(x) || !is.numeric(x))
	stop("'x' must be 2d numeric matrix")
	if (scale) x <- scale(x)
	 d <- dist(x,method=dmethod)
	if (missing(ser_weight))
	 if (isTRUE(all.equal(cost, costLS)))
	    ser_weight <- x[,1]
    	 else ser_weight <- as.matrix(d)

     dser.dist(d,ser_weight,cost=cost,...)
	}
	
dser.dist <- function(x,ser_weight,cost=costBAR,hmethod="average",...) {
    h <- hclust(x,method=hmethod)
    if (missing(ser_weight) && !isTRUE(all.equal(cost, costLS)))
    	ser_weight <- as.matrix(x)
    DendSer(h,ser_weight,cost=cost,...)
	}
	
dser.hclust <- function(x,ser_weight,cost=costBAR,...) {
     DendSer(x,ser_weight,cost=cost,...)
	}

