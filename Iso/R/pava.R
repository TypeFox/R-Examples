pava <- function(y, w=NULL, decreasing=FALSE, long.out = FALSE, stepfun=FALSE)
{
#
# Function 'pava'.  To perform isotonic regression for a simple
# (increasing) linear ordering using the ``pool adjacent violators
# algorithm''.  If long.out = TRUE then the result returned consists
# of a list containing the fitted values, the final weights, and a set
# of indices `tr', made up of the smallest index in each level set,
# which thus keeps track of the level sets.  Otherwise only the fitted
# values are returned.
# 
	if(decreasing) y <- rev(y)
	n <- length(y)
	if(is.null(w))
		w <- rep(1, n)
	else if(decreasing) w <- rev(w)
	if(n == 1) {
		if(long.out) return(list(y=y,w=w,tr=1))
		else return(y)
	}

	rslt <- .Fortran(
			"pava",
			y=as.double(y),
			w=as.double(w),
			kt=integer(n),
			n=as.integer(n),
			PACKAGE="Iso"
			)
	y <- if(decreasing) rev(rslt$y) else rslt$y
	if(long.out | stepfun ) {
		tr <- rslt$kt
		if(decreasing) 
			tr <- unname(unlist(tapply(1:n,-rev(tr),
				function(x){rep(min(x),length(x))})))
	}
	if(long.out) {
		w <- if(decreasing) rev(rslt$w) else rslt$w
		lout <- list(y = y, w = w, tr = tr)
	}
	if(stepfun) {
		knots <- 1+which(diff(tr)!=0)
		y0    <- c(y[1],y[knots])
		h     <- stepfun(knots,y0)
	}
	ntype <- 1+sum(c(long.out,stepfun)*(1:2))
	switch(ntype,y,lout,h,c(lout,list(h=h)))
}
