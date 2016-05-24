"CqN" <- 
function(x, q=2, method="pi.a"){
	if(q < 2){ stop("The function is currently not defined at q < 2") }
	N <- nrow(x)
	# check out the version:
	x.pi <- switch(method, pi.a = x/rowSums(x), pi.g = colSums(x/sum(x)), pi.ga = x/sum(x))
	# calculate the numerator (see Chao et al. 2008)
	num <- sum(apply(x.pi, 2, function(x) sum(x, na.rm=TRUE)^q - sum(x^q, na.rm=TRUE)))/(N^q-N)
	# calculate the denominator (see Chao et al. 2008)
	den <- sum(apply(x.pi, 2, function(x) sum(x^q, na.rm=TRUE)))/N
	CqN <- num/den
	names(CqN) <- paste("C", q, N, sep="")
	return(CqN)
}