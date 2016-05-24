
########################################################
# extract replicated parameters for BIFIE method
extract.replicated.pars <- function( BIFIE.method , type=NULL ){
    
	parsM <- parsrepM <- NULL
	res1 <- BIFIE.method
	#************************************
	#**** linear regression
	if ( class( BIFIE.method) == "BIFIE.linreg"){ 	
		# parameters in every imputed dataset
		parsM <- res1$output$regrcoefM
		# replicated parameters
		parsrepM <- res1$output$regrcoefrepM
				}
	#************************************
	#**** path model
	if ( class( BIFIE.method) == "BIFIE.pathmodel"){ 	
		# parameters in every imputed dataset
		parsM <- res1$output$parsM
		# replicated parameters
		parsrepM <- res1$output$parsrepM
				}				
	#************************************
	#**** correlation
	if ( class( BIFIE.method) == "BIFIE.correl"){ 	
		parsM <- res1$output$cor1M
		parsrepM <- res1$output$cor1repM
		if ( ! is.null(type) ){
			if ( type == "cov"){		
				parsM <- res1$output$cov1M
				parsrepM <- res1$output$cov1repM		
								}
							}
				}				
	#************************************
	#**** frequencies
	if ( class( BIFIE.method) == "BIFIE.freq"){ 	
		parsM <- res1$output$perc2M
		parsrepM <- res1$output$perc2repM
				}		
	#************************************
	#**** univar
	if ( class( BIFIE.method) == "BIFIE.univar"){ 	
		parsM <- res1$output$mean1M
		parsrepM <- res1$output$mean1repM
				}		
	#************************************
	#**** crosstab
	if ( class( BIFIE.method) == "BIFIE.crosstab"){ 	
		parsM <- res1$output$ctparsM
		parsrepM <- res1$output$ctparsrepM
				}
	#************************************
	#**** logistreg
	if ( class( BIFIE.method) == "BIFIE.logistreg"){ 	
		parsM <- res1$output$regrcoefM
		parsrepM <- res1$output$regrcoefrepM
				}							
	#************************************
	#**** BIFIE.by
	if ( class( BIFIE.method) == "BIFIE.by"){ 	
		parsM <- res1$output$parsM
		parsrepM <- res1$output$parsrepM
				}		
	#************************************
	#**** BIFIE.derivedParameters
	if ( class( BIFIE.method) == "BIFIE.derivedParameters"){ 
		parsM <- res1$parsM
		parsrepM <- res1$parsrepM
				}						
	#************************************
	#**** BIFIE.twolevelreg
	if ( class( BIFIE.method) == "BIFIE.twolevelreg"){ 	
		parsM <- res1$output$parsM
		parsrepM <- res1$output$parsrepM
				}						
	res <- list( "parsM" = parsM , "parsrepM" = parsrepM ,
		"parnames"= res1$parnames )
	return(res)
		}
					