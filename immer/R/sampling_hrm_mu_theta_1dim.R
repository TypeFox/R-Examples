
##########################################################
# sampling theta one dimension
sampling_hrm_mu_1dim <- function( theta , prior , N ){
	m1 <- mean( theta )
	m2 <- prior$mu$M
	w1 <- N / stats::var(theta) 
	w2 <- 1 / prior$mu$SD^2
	m0 <- ( w1*m1 + w2*m2 ) / (w1 + w2 )
	s0 <- 1 / sqrt( w1 + w2 )
	mu_new <- stats::rnorm( 1 , mean = m0 , sd = s0 )
	return(mu_new)
		}
###########################################################		