
######################################################################
IRT.derivedParameters <- function( jkobject , derived.parameters ){
		object <- jkobject
		parsM <- object$parsM
		est <- object$jpartable$value
		names(est) <- rownames(parsM)
		parsM <- t(parsM)
		allformulas <- derived.parameters[[1]]
		FF <- length(derived.parameters)
		if (FF>1){
		for (ff in 2:FF){
			t1 <- stats::terms( allformulas)    
			t2 <- paste( c( attr( t1 , "term.labels" ) , 
				attr(  stats::terms( derived.parameters[[ff]] ) , "term.labels" )  ),
				collapse= " + " )
			allformulas  <- stats::as.formula( paste( " ~ 0 + " , t2 ) )
					}
					}
		# create matrices of derived parameters
		der.pars <- stats::model.matrix( allformulas , as.data.frame( t(est) ) )
		colnames(der.pars) <- names(derived.parameters)
		der.pars.rep <- stats::model.matrix( allformulas , as.data.frame( parsM) )
		colnames(der.pars.rep) <- names(derived.parameters)
		parsM <- t(der.pars.rep)
		res0 <- jkestimates( est = as.vector(der.pars) , parsM = parsM , fayfac = object$fayfac )
		jpartable <- data.frame( "parnames" = names(derived.parameters) , 
						"value" = as.vector(der.pars) )
		jpartable$jkest <- res0$dfr$jkest		
		jpartable$jkse <- res0$dfr$jkse									
		res <- list( "parsM" = parsM , "vcov" = res0$vcov_pars , 
				      "jpartable" = jpartable , "fayfac" = object$fayfac )
		class(res) <- "IRT.jackknife"
		return(res)					
				}	
###############################################################