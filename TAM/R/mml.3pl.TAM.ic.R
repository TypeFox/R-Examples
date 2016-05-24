

##################################
# Information criteria
.mml.3pl.TAM.ic <- function( nstud , deviance , xsi , xsi.fixed ,
	beta , beta.fixed , ndim , variance.fixed , G , irtmodel ,
	B_orig=NULL , B.fixed , E , est.variance , resp ,
	est.slopegroups=NULL , skillspace , delta , delta.fixed , est.guess , fulldesign ,
	est.some.slopes , gammaslope , gammaslope.fixed, gammaslope.constr.V ,
	gammaslope.constr.Npars , gammaslope.center.index ,
    gammaslope.prior , numdiff.parm ){
  #***Model parameters
  h <- numdiff.parm
  ic <- data.frame("n" = nstud , "deviance" = deviance )
  dev <- deviance
	# xsi parameters
	ic$Nparsxsi <- length(xsi)
	if ( ! is.null( xsi.fixed) ){ 
			ic$Nparsxsi <- ic$Nparsxsi - nrow(xsi.fixed ) }
	# B slopes
	ic$NparsB <- 0
    if ( est.some.slopes ){
		ic$NparsB <- length(gammaslope)
		if ( ! is.null(gammaslope.constr.V) ){
		   ic$NparsB <- ic$NparsB - ncol(gammaslope.constr.V)		
							}
						}
	if ( ! is.null( gammaslope.fixed ) ){
		ic$NparsB <- ic$NparsB - nrow(gammaslope.fixed )
					}
	ic$NparsB <- ic$NparsB - gammaslope.constr.Npars
	if ( ! is.null(gammaslope.center.index ) ){
		ic$NparsB <- ic$NparsB - max( gammaslope.center.index )
						}
    # non-active gammaslope parameters
	ic$Ngamma.nonactive <- 0
	if ( ! is.null(gammaslope.prior ) ){
	     if ( ncol(gammaslope.prior) > 2 ){
			  ic$Ngamma.nonactive <- ic$Ngamma.nonactive + 
						sum( gammaslope < gammaslope.prior[,3] + 3*h )
			  ic$Ngamma.nonactive <- ic$Ngamma.nonactive + 
						sum( gammaslope > gammaslope.prior[,4] - 3*h )	
							}
					}
						
	# beta regression parameters
	ic$Nparsbeta <- 0
	ic$Nparscov <- 0
	if ( skillspace=="normal"){
		ic$Nparsbeta <- dim(beta)[1] * dim(beta)[2]
		if ( ! is.null( beta.fixed) ){ 
				ic$Nparsbeta <- ic$Nparsbeta - nrow(beta.fixed ) }
		# variance/covariance matrix
		ic$Nparscov <- ndim + ndim*(ndim-1)/2
		if ( ! est.variance ){ ic$Nparscov <- ic$Nparscov - ndim }
		if ( ! is.null( variance.fixed) ){ 
				ic$Nparscov <- max(0 , ic$Nparscov - nrow(variance.fixed ) )
										 }
								}												 
	ic$Nguess <- length( setdiff( unique(est.guess) , 0 ) )									 
	if ( skillspace != "normal" ){								 
		ic$Ndelta <- prod( dim(delta) )
#        if ( fulldesign ){ 
             ic$Ndelta <- ic$Ndelta - ncol(delta)
#							}
             if ( ! is.null( delta.fixed ) ){							
				ic$Ndelta <- ic$Ndelta - nrow(delta.fixed )
									}			 							
							} else { 
			ic$Ndelta <- 0 
						}
	# total number of parameters
	ic$Npars <- ic$np <- ic$Nparsxsi + ic$NparsB + ic$Nparsbeta + ic$Nparscov + 
		          ic$Nguess + ic$Ndelta	- ic$Ngamma.nonactive
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