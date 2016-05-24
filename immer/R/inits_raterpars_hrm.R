
###########################################
# inits rater parameters
inits_raterpars_hrm <- function( rater , I , est_settings ){
	R <- length( unique(rater) )
	phi <- matrix( 0 , nrow=I , ncol=R )
	psi <- .3 + phi
	if ( est_settings$est.psi == "n" ){
		psi <- 1E-10 + 0*phi
					}
	res <- list(R=R , phi=phi , psi=psi)
	return(res)
			}
#############################################			