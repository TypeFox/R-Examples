
#################################################################
# categorize variables into classes
categorize <- function( dat , categorical=NULL , quant=NULL ,
			lowest=0){

	dat2 <- dat

	#*********************
	# categorize variables, numbered from lowest to max
	dfr <- NULL
	if (! is.null( categorical) ){
		VV <- length(categorical)	
		for (vv in 1:VV){
			# vv <- 1
			var.vv <- categorical[vv]    
			dat.vv <- dat[,var.vv]
			vals.vv <- base::sort( unique( dat.vv ) )
			dfr.vv <- data.frame( "index" = vv , "variable" = var.vv ,
					"column" = which( colnames(dat) == var.vv ) ,
					"orig" = vals.vv , 
					"recode" = seq( lowest , length(vals.vv) -1 + lowest)
								)      
			dfr <- base::rbind( dfr, dfr.vv )
			dat2[ , var.vv] <- base::match( dat[,var.vv] , vals.vv ) - 1 + lowest
					}
					}
	#**************************
	# categorize continuous variables into quantiles
	
	dfr2 <- NULL
	
	if ( ! is.null(quant) ){
		vars <- names(quant)
		VV <- length(vars)
	
		for (vv in 1:VV){
			# vv <- 1
			vars.vv <- vars[vv]
			q1 <- quant[ vars.vv ]
			quant.vv <- stats::quantile( dat[,vars.vv] , na.rm=TRUE , 
						   prob = seq( 0 , 1 , len=q1+1 )  )
			quant.vv[1] <- quant.vv[1] - 1
			quant.vv[q1+1] <- quant.vv[q1+1] + 1		
			quant.vv <- base::unique( quant.vv )

			m1 <- base::cut( dat[,vars.vv]  , breaks = quant.vv )
			m2 <- base::sort( unique(m1) )

			dfr2.vv <- data.frame( "index" = vv , "variable" = vars.vv ,
						"column" = which( colnames(dat) == vars.vv ) ,
						"orig" = m2 , 
						"recode" = seq( 0 , length(m2) -1 + lowest)
									)  				
			dfr2 <- base::rbind( dfr2 , dfr2.vv )							
			dat2[ , vars.vv ] <- base::match( m1 , m2 ) - 1 + lowest		
				}

						}
						
	#***************************
	# OUTPUT
	res <- list( "data" = dat2 , "categ_design" = dfr ,
			"quant_design" = dfr2 )
	return(res)	
}
#################################################################			
