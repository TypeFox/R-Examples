"dnormt3" <-
function(x, mean=0, sd=1){

if (sd==0)	{
	f <- rep(0, length(x))
	f[x==0] <- Inf
	}
else	{
	f <- sqrt(2)*normt3ip(mu=(x-mean)*sqrt(2)/sd, sigma=1)/sd
	sv <- abs(x - mean)*sqrt(2)/sd >= 37
	f[sv] <- 0
	}
f
}
