ed1 <-
function(x, q=1, retq = TRUE) {
	bsums <- function(x, q=1) {
		if (all(x==0)) return(0)
		if (q==0) sum(x>0)
		else if (q==1) sum(-x*log(x),na.rm=TRUE)
		else sum(x^q)
	}
	if (!is.null(dim(x))) a <- apply(x,1,ed1,q=q,retq=retq)
	else {
		x <- matrix(x,nrow=1)
		rs <- sum(x)
		if (rs==0) stop("Zero data")
		x <- x/rs
		a <- bsums(x,q=q)
		if (retq) {
			if (q==1) {
				a <- exp(a)
			}
			else if (q!=1) {
				a <- a^(1/(1-q))
			}
		}
	}
	a
}

