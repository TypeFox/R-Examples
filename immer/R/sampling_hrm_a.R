
#####################################################################
# sampling of item parameters
sampling_hrm_a <- function( xi , xi_ind  , b , a , maxK , prior , MHprop , I , 
		theta , eps=1E-20 ){
	# refresh count
	MHprop$refresh_count$a <- MHprop$refresh_count$a + 1	
	for (ii in 1:I ){
			a_new <- a_old <- a
			a_new[ii] <- stats::rlnorm( 1 , meanlog= log(a_old[ii]) , sdlog = MHprop$SD$a[ii] )			
			p_new <- stats::dlnorm( a_new[ii] , meanlog = prior$a$M[ii] , sdlog = prior$a$SD[ii] ) 
			p_old <- dlnorm( a_old[ii] , meanlog = prior$a$M[ii] , sdlog = prior$a$SD[ii] ) 				
			ll_new <- sum( log( probs_gpcm( x = xi[,ii] , theta , b=b[ii,] , 
			         a=a_new[ii] , K=maxK[ii] , x_ind = xi_ind[,ii] )  + eps) )
			ll_old  <- sum( log( probs_gpcm( x = xi[,ii] , theta , b=b[ii,] , 
			          a=a_old[ii] , K=maxK[ii] , x_ind = xi_ind[,ii] )  + eps) )
			ratio <- p_new * exp( ll_new - ll_old ) / ( p_old  )
			if ( ratio > stats::runif(1) ){
				MHprop$accept$a[ii] <- MHprop$accept$a[ii] + 1 
				a <- a_new
								 }
						}  # end item ii
        res <- list( a = a , MHprop = MHprop )
        return(res)
				}
#####################################################################