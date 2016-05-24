
#####################################
# sampling from a truncated normal distribution
rtrnorm <- function( N , mean , sd , lower = rep(-Inf,N) , upper=rep(Inf,N) ){
	t1 <- stats::pnorm( lower , mean = mean , sd = sd )
	t2 <- stats::pnorm( upper , mean = mean , sd = sd )
	rn <- stats::runif( N , t1 , t2 )
	stats::qnorm( rn , mean = mean , sd = sd )
				}
#########################################				