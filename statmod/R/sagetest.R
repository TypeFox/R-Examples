#  SAGE.R

sage.test <- function(x, y, n1=sum(x), n2=sum(y))
#	Exact binomial probabilities for comparing SAGE libraries
#	Gordon Smyth
#	15 Nov 2003.  Last modified 20 July 2012.
{
	if(any(is.na(x)) || any(is.na(y))) stop("missing values not allowed")
	x <- round(x)
	y <- round(y)
	if(any(x<0) || any(y<0)) stop("x and y must be non-negative")
	if(length(x) != length(y)) stop("x and y must have same length")
	n1 <- round(n1)
	n2 <- round(n2)
	if(!missing(n1) && any(x>n1)) stop("x cannot be greater than n1")
	if(!missing(n2) && any(y>n2)) stop("y cannot be greater than n2")
	size <- x+y
	p.value <- rep(1,length(x))
	if(n1==n2) {
		i <- (size>0)
		if(any(i)) {
			x <- pmin(x[i],y[i])
			size <- size[i]
			p.value[i] <- pmin(2*pbinom(x,size=size,prob=0.5),1)
		}
		return(p.value)
	}
	prob <- n1/(n1+n2)
	if(any(big <- size>10000)) {
		ibig <- (1:length(x))[big]
		for (i in ibig) p.value[i] <- chisq.test(matrix(c(x[i],y[i],n1-x[i],n2-y[i]),2,2))$p.value
	}
	size0 <- size[size>0 & !big]
	if(length(size0)) for (isize in unique(size0)) {
		i <- (size==isize)
		p <- dbinom(0:isize,prob=prob,size=isize)
		o <- order(p)
		cumsump <- cumsum(p[o])[order(o)]
		p.value[i] <- cumsump[x[i]+1]
	}
	p.value
}
