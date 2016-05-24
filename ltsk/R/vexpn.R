vexpn <-
function(x,sill,range,nugget)
{
## note 07/13/2015
## this is Okay as long as x > 0

  if( range==0) return(rep(sill+nugget,length(x)))
  nugget + sill *( 1 - exp(-3*x/range))
}
cexpn <- function(x,sill,range,nugget)
{
## this allows x to be small relative to range
	totvar <- sill+nugget ## total variance
	r <- rep(totvar,length(x)) ## return vector
	chk <- (abs(x)/range) > 1e-7 ## non-zero lag
	lx <- x[chk] ## 
	lc <- sill * exp(-3*lx/range) ## covariance
	r[chk] <- lc ## storage
	r ## return variance
}
