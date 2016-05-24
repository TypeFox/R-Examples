cca.single <- function(data, s, x,y, mode=3){
	#do checks
	stopifnot(is.numeric(data))
	stopifnot(is.matrix(data))
	stopifnot(is.numeric(s))
	stopifnot(is.numeric(x))
	stopifnot(is.numeric(y))
	x <- as.integer(x) - 1 # as C starts counting at 0
	y <- as.integer(y) - 1 # as C starts counting at 0
	xmax <- nrow(data) 
	ymax <- ncol(data) 
	stopifnot(x>=0 & x<xmax) # starting cell must be inside
	stopifnot(y>=0 & y<ymax) # starting cell must be inside
	stopifnot(mode==1 | mode==2 | mode==3)
	method <- switch(mode, "burnn", "burns", "burnr")
	the.data <- as.integer(t(data))
	clu <- as.integer(rep(0, ncol(data)*nrow(data)))
	count.max <- as.integer(ncol(data)*3)
	count <- as.integer(rep(0, count.max))
	if(mode==1){
		out <- .C(method, x=as.integer(x), y=as.integer(y),  c=as.integer(1), xmax=as.integer(xmax), ymax=as.integer(ymax),  data=the.data, clu=clu, CLASSES=c("integer", "integer", "integer", "integer", "integer", "integer","integer"))
	} else {
		out <- .C(method,  data=as.integer(the.data), clu=as.integer(clu), x=as.integer(x), y=as.integer(y),  c=as.integer(1), s=as.integer(s), xmax=as.integer(xmax), ymax=as.integer(ymax))
	}
	
	return(matrix(out$clu, ncol=ncol(data), byrow=TRUE))
}


