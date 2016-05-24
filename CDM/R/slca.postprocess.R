


#############################################################
# calculation of information criteria and number of parameters
.slca.calc.ic <- function( dev , dat , G ,   
			K, TP ,I , delta.designmatrix , delta.fixed ,
			Xlambda , Xlambda.fixed , data0 , deltaNULL ,
			Xlambda.constr.V 
				){
    ic <- list( "deviance" = dev , "n" = nrow(data0) )
	ic$traitpars <- 0
	ic$itempars <- 0	
	
	ic$itempars <- length(Xlambda)
	if ( ! is.null(Xlambda.fixed ) ){
		ic$itempars <- ic$itempars - nrow(Xlambda.fixed )
									}
	if ( ! is.null( Xlambda.constr.V ) ){
		ic$itempars <- ic$itempars - ncol(Xlambda.constr.V )
									}
																
	ic$traitpars <- G * ncol(delta.designmatrix ) - G*deltaNULL
	if ( ! is.null(delta.fixed ) ){
		ic$traitpars <- ic$traitpars - nrow(delta.fixed )
									}
	#***********************************************
	# information criteria
	ic$np <- ic$itempars + ic$traitpars	
#	ic$n <- n # number of persons
	# AIC
        ic$AIC <- dev + 2*ic$np
        # BIC
        ic$BIC <- dev + ( log(ic$n) )*ic$np
        # CAIC (conistent AIC)
        ic$CAIC <- dev + ( log(ic$n) + 1 )*ic$np
		# corrected AIC
        ic$AICc <- ic$AIC + 2*ic$np * ( ic$np + 1 ) / ( ic$n - ic$np - 1 )				
    return(ic)
		}
###################################################################
		