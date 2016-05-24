
##############################################################
lsem.fitsem <- function( dat , weights , lavfit ,
			fit_measures , NF , G , moderator.grid , verbose ,
			pars , standardized ){

    parameters <- NULL
	fits <- NULL 	
	pars0 <- pars	
	
	if (verbose){
		cat( "** Fit lavaan model\n")
		G1 <- min(G,10)	
		pr <- round( seq(1,G , len=G1) )
		cat("|")
		cat( paste0( rep("*",G1) , collapse="") )
		cat("|\n")
		cat("|")
				}
			

	for (gg in 1:G){
		# gg <- 1			
		dat$weight <- weights[,gg]
		datsvy <- survey::svydesign(id=~index,   weights=~weight ,    data=dat)
		# fit the model using weighted data
		survey.fit <- lavaan.survey::lavaan.survey(lavfit, datsvy )
		dfr.gg <- pars <- lavaan::parameterEstimates(survey.fit) 		
		if (standardized){			
			sol <- lavaan::standardizedSolution( survey.fit )
			colnames(sol)[ which( colnames(sol) == "est.std" ) ] <- "est"
			sol$lhs <- paste0( "std__" , sol$lhs)
			pars <- plyr::rbind.fill( pars , sol )	
			dfr.gg <- pars
						} 
						
		pars <- paste0( pars$lhs , pars$op , pars$rhs )					
		NP <- length(pars0)
		ind <- match( pars0 , pars )
		dfr.gg <- dfr.gg[ ind , ]
		dfr.gg <- data.frame("grid_index"=gg , "moderator" = moderator.grid[gg] ,
						  "par"= pars0 , "parindex" = 1:NP , dfr.gg	)
		dfr.gg0 <- data.frame("grid_index"=gg , "moderator" = moderator.grid[gg] ,
						  "par"= fit_measures , "parindex" = NP + 1:NF , 
						  "est"= lavaan::fitMeasures(survey.fit , fit.measures= fit_measures ) ,
						  "op"="fit" )
		vars <- setdiff( colnames(dfr.gg) , colnames(dfr.gg0) )
		for (vv in vars){ dfr.gg0[,vv] <- NA }
		dfr.gg <- rbind( dfr.gg , dfr.gg0[ , colnames(dfr.gg) ] )		
		parameters <- rbind( parameters , dfr.gg ) 
		# fits <- rbind( fits , dfr.gg ) 
		if (verbose){
			if ( gg %in% pr ){
				cat("-") ; utils::flush.console()
					}
					}
				}
	if (verbose){
		cat("|\n")
		utils::flush.console()
				}


	parameters <- parameters[ order(parameters$parindex) , ]	
#	fits <- fits[ order(fits$fitindex) , ]	
#	rownames(fits) <- NULL
				
	res <- list( "parameters" = parameters ) #  , "fits" = fits )
	return(res)
				
			}
#######################################################################			