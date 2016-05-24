

#################################################################
# approximation of the bivariate normal integral
# by the method of Cox and Wermuth (1991)
normal2.cw <- function( a , b , rho ){
    # (X,Y) ~ N_2(0 , rho)
    #   P( X > a1 , Y > b1 , rho ) 
    # = P( -X < - a1 , -Y < -b1 , rho )
	if ( any( rho < 0 )){stop("only positive correlations are allowed!")}
    a11 <- a1 <- - a
    b1 <- - b
#	# APPROXIMATION for negative correlations
#	ind.neg <- which( rho < 0 )
#	if ( length(ind.neg) > 0 ){
#		rho[ind.neg] <- - rho[ind.neg]
#		b1[ind.neg] <- -b1[ind.neg]
#					}
	# APPROX. (ii)
	ind2 <- which( a1 < 0 & b1 < 0 )
	if ( length(ind2) > 0 ){
		a1[ ind2 ] <- - a1[ind2] 
		b1[ ind2 ] <- - b1[ind2] 
				}
	# APPROX. (i):
	# a > 0 and b > 0 => a > b
	ind <- which( a1 < b1 )	
	t1 <- b1
	if ( length(ind) > 0 ){ 
		b1[ ind ] <- a1[ind]
		a1[ind] <- t1[ind]
			}
    t1 <- stats::pnorm( - a1 )
    mu <- stats::dnorm( a1 ) / stats::pnorm( - a1 )
    xi <-  ( rho * mu - b1 ) / sqrt( 1 - rho^2 ) 
    # see Hong (1999)
    # sig2 <- 1 + mu - mu^2 # => formula in Cox & Wermuth (1991)
    sig2 <- 1 + a1*mu - mu^2
    #  Formula (3)
    # t1 * pnorm( xi )
    # Formula (4): correction in Hong (1999)
    prob1 <- t1 * ( stats::pnorm(  xi ) - 1/2 * rho^2 / ( 1 - rho^2 ) * xi * 
				stats::dnorm( xi ) * sig2  )
	# adjust formula in case of APPROX. (ii)
	if ( length(ind2) > 0 ){
		# CW. Formula in (ii), p. 264
		prob1[ind2] <- 1 - stats::pnorm( -a1[ind2] ) - stats::pnorm( -b1[ind2] ) + prob1[ind2] 
			}
#	# negative correlations
#	if ( length(ind.neg) > 0 ){
#		prob1[ind.neg] <- pnorm( - a11[ind.neg] ) - prob1[ ind.neg]
#					}	
    return( prob1 )
        }
#################################################################

