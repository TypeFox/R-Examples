weibull.int <-
function(v.dat, msg) {
### internal function for calculation weibull parameters
 		
 	if(any(is.na(v.dat)==TRUE)) {
 		n <- length(v.dat[is.na(v.dat)==TRUE])
 		v.dat<- v.dat[!is.na(v.dat)]
 		if(msg) message(n, " NA found and excluded from calculation")
 	}
	if(any(v.dat<=0)) {
		n <- length(v.dat[v.dat<=0])
 		v.dat<- v.dat[!is.na(v.dat)]
		if(msg) message(n, " none-positives found and excluded from calculation")
	}
	m1 <- mean(v.dat)
	m3 <- mean(v.dat^3)
	ratio.data <- m1^3 / m3
	
	k.theo <- seq(0.02, 10, 0.01)
	ratio.theo <- gamma(1+1/k.theo)^3 / gamma(1+3/k.theo) 	
	ind.k <- max(which(ratio.theo<ratio.data))
	
	delta.k <- k.theo[ind.k+1] - k.theo[ind.k]
	delta.ratio <- ratio.theo[ind.k+1] - ratio.theo[ind.k]
	delta.mess  <- ratio.data - ratio.theo[ind.k]
	
	k <- k.theo[ind.k] + delta.k / delta.ratio * delta.mess
	A <- m1 / gamma(1+1/k)
		
	return(list(A=A, k=k))
}
