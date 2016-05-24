
########################################
# latent regression information criteria
latreg_TAM.ic <- function( nstud , deviance , 
	beta , beta.fixed , ndim , variance.fixed , G , 
	est.variance , variance.Npars=NULL , group ){

  #***Model parameters
  ic <- data.frame("n" = nstud , "deviance" = deviance )
  dev <- deviance
	# xsi parameters
	ic$Nparsxsi <- 0
	# B slopes
	ic$NparsB <- 0	
	# beta regression parameters
	ic$Nparsbeta <- dim(beta)[1] * dim(beta)[2]
	if ( ! is.null( beta.fixed) ){ 
			ic$Nparsbeta <- ic$Nparsbeta - nrow(beta.fixed ) }
	# variance/covariance matrix
	ic$Nparscov <- ndim + ndim*(ndim-1)/2
	if ( ! est.variance ){ ic$Nparscov <- ic$Nparscov - ndim }
	if ( ! is.null( variance.fixed) ){ 
			ic$Nparscov <- max(0 , ic$Nparscov - nrow(variance.fixed ) )
									 }
	if ( ! is.null(variance.Npars) ){
	       ic$Nparscov <- variance.Npars 
						}
	if ( ! is.null(group) ){
	     ic$Nparscov <- ic$Nparscov + length( unique(group) ) - 1
							}
	# total number of parameters
	ic$Npars <- ic$np <- ic$Nparsxsi + ic$NparsB + ic$Nparsbeta + ic$Nparscov
    	# AIC
        ic$AIC <- dev + 2*ic$np
		# AIC3
		ic$AIC <- dev + 3*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
		# adjusted BIC 
		ic$aBIC <- dev + ( log( ( ic$n -2 ) / 24 ) )*ic$np
        # CAIC (consistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )	  
		
	return(ic)
	}