valfpa <- function(vw,b,delta,funcpa){
	bk <- b + (exp(vw)*delta)
	return(-funcpa(bk))
}