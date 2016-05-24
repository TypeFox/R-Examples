
######################################################
# calculate degrees of freedom according to Rubin
rubin_calc_df <- function( res_pars , Nimp , indices=NULL , digits=2){
    W <- res_pars$pars_varWithin
    B <- res_pars$pars_varBetween	
	df <- rubin_calc_df2( B , W , Nimp , indices=indices , digits=digits)	
    return(df)
        }
######################################################		


######################################################
rubin_calc_df2 <- function( B , W , Nimp , indices=NULL , digits=2){
	if ( ! is.null(indices) ){
		W <- W[indices]
		B <- B[indices]
			}
	B <- B + W * 1E-15
    df <- ( 1 + Nimp / ( Nimp + 1 ) * W / B  )^2 * (Nimp - 1 )
	df <- round( df , digits )
	df <- ifelse( df > 1000 , Inf , df )
    return(df)
        }
######################################################	