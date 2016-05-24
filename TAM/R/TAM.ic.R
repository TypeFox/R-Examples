

##################################
# Information criteria
.TAM.ic <- function( nstud , deviance , xsi , xsi.fixed ,
	beta , beta.fixed , ndim , variance.fixed , G , irtmodel ,
	B_orig=NULL , B.fixed , E , est.variance , resp ,
		est.slopegroups=NULL , variance.Npars=NULL , group ){
	# 2PL estimation
	# c("2PL","GPCM","GPCM.design","2PL.groups") )	
# Revalpr( "unique(group)")
  #***Model parameters
#  nstud <- sum( 1*(rowSums( 1 - is.na(resp) ) > 0  ) )
  ic <- data.frame("n" = nstud , "deviance" = deviance )
  dev <- deviance
	# xsi parameters
	ic$Nparsxsi <- length(xsi)
	if ( ! is.null( xsi.fixed) ){ 
			ic$Nparsxsi <- ic$Nparsxsi - nrow(xsi.fixed ) }
	# B slopes
	ic$NparsB <- 0
	if ( irtmodel == "2PL" ){
		ic$NparsB <- sum( B_orig != 0 )
						    }
	if ( irtmodel == "GPCM" ){
		ic$NparsB <- ncol(resp)
						    }
	if ( irtmodel == "GPCM.design" ){
		ic$NparsB <- ncol(E)
						    }							
	if ( irtmodel == "2PL.groups" ){
#		ic$NparsB <- sum( B_orig != 0 )
		ic$NparsB <- length( unique( est.slopegroups ) )
		# This is not yet correct for multiple dimensions and multiple
		# categories
						    }								
	if ( ! is.null( B.fixed ) ){
		ic$NparsB <- ic$NparsB - nrow(B.fixed )
					}
	
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
		ic$AIC3 <- dev + 3*ic$np
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