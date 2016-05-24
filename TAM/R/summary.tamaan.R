
###############################################
# summary tamaan
summary.tamaan <- function( object , file=NULL , ... ){
	if ( ! is.null( file ) ){
			sink( paste0( file , "__SUMMARY.Rout") , split=TRUE )
						}	
	#**********************
	# general tamaan syntax
	# cat("------------------------------------------------------------\n")
	cat("!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:!:\n")
	cat(paste0("tamaan function using '" ,
			object$tamaanify$method , "' method\n\n") )
	cat( paste(object$tamaanify$tammodel) )
	cat("\n\n")

	#**********************
	# tam.mml
	if ( object$tamaan.method == "tam.mml" ){
		summary.tam.mml( object , file=NULL ,... )
							}
	#**********************
	# tam.mml.2pl
	if ( object$tamaan.method == "tam.mml.2pl" ){
		summary.tam.mml( object , file=NULL , ... )
							}
    #**********************
	# tam.mml.3pl
	if ( object$tamaan.method == "tam.mml.3pl" ){
		# overview for all parameters
		summary.tamaan.3pl.intro(object)
		
		# distribution discrete skill space
		if ( !( object$tamaanify$ANALYSIS.list$type %in% c( "MIXTURE" ) )){
			if ( object$skillspace == "discrete" ){	
				summary.tamaan.3pl.discrete.distribution(object)		
								}
					}
		# distribution normal skill space
		if (object$skillspace == "normal"){
			summary.tamaan.normal.skillspace(object)
									}
			
		# cluster locations		
		if ( object$tamaanify$ANALYSIS.list$type %in% c( "LOCLCA" ) ){
			summary.tamaan.3pl.loclca(object)
								}

		# distribution mixture	
		if ( object$tamaanify$ANALYSIS.list$type %in% c( "MIXTURE" ) ){
			summary.tamaan.3pl.distr.mixture(object)
								}
								
		# Item parameters
		print_ipars <- FALSE
		if (object$skillspace == "normal"){ print_ipars <- TRUE }
		if ( object$tamaanify$ANALYSIS.list$type %in% c( "LOCLCA" ) ){
				print_ipars <- TRUE 
							}
				
		
		if ( print_ipars ){
			summary.tamaan.item.parameters(object)
									}
				
		# Latent class probabilities
		if ( object$tamaanify$ANALYSIS.list$type %in% c( "LCA" , "OLCA" ) ){
			summary.tamaan.3pl.lcaprobs(object)
								}
								
		# Mixture distribution
		if ( object$tamaanify$ANALYSIS.list$type %in% c( "MIXTURE" ) ){
				summary.tamaan.item.parameters.mixture( object )
							}
		
								
		
					}
	#**********************************	
		if ( ! is.null( file ) ){
			sink(  )
							}
				}
################################################				
