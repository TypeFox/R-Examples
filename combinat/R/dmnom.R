"dmnom"<-
function(x, size = sum(x), prob = stop("no prob arg"))
{
#       DATE WRITTEN: 22 May 1995           LAST REVISED:  22 May 1995
#       AUTHOR:  Scott Chasalow
#
	p <- max(length(x), length(prob))
	x <- rep(x, length = p)
	prob <- rep(prob, length = p)
	prob <- prob/sum(prob)
	if(sum(x) != size)
		0
	else exp(logfact(size) + sum(x * log(prob) - logfact(x)))
}

