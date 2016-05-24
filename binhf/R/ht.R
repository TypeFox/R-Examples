"ht" <-
function (x) 
{
	n <- length(x)
	i <- floor(log(n, base=2))
	d <- c();s<- c()
	#detail
	for (j in 1:i){
		x1 <- as.vector(na.omit(filter(x, c(rep(-1/{2^j}, 2^{j-1}), rep(1/{2^j}, 2^{j-1})) )))
		x2 <- as.vector(na.omit(filter(x, c(rep(1/{2^j}, 2^{j-1}), rep(1/{2^j}, 2^{j-1})) )))
		x1 <- x1[seq(1,length(x1), by=2^j)]
		x2 <- x2[seq(1,length(x2), by=2^j)]
		s <- c(s,x2)
		d <- c(d,x1)
	}
	#smooooth - just s_0
#	s <- as.vector(na.omit(filter(x, c(1/2, 1/2) )))
#	s <- s[seq(1,length(s), by=2)]
	hf <- list(d=d,s=s)
	return(hf)
}

