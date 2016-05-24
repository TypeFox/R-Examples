
############################################################
mfr.dataprep <- function( formulaA , xsi.setnull , B , Q ,
		resp, pid, facets , beta.fixed ){
	
	tA <- stats::terms( formulaA )
	tlab <- attr(tA , "factors")
	
	# redefine formula
	stlab <- apply( tlab , 1 , sum )		
	ind <- which( stlab == 0 )
	nullfacets <- names(stlab)[ind]

	if ( is.null(pid) ){
		pid <- seq( 1 , nrow(resp) )
		cat("--- Created person identifiers.\n")
					}
		

	#********************
	# data restructuring for non-identifiable combinations
	facets_labs <- setdiff( rownames(tlab) , c("item" , "step") )
	
	# create combination of pid and facets
	FF <- length(facets_labs)	
	combi <- pid
	if ( FF>0){
		for (ff in 1:FF){
			combi <- paste0( combi , "-" , facets[ , facets_labs[ff] ] )
						}
					}

	dups <- base::duplicated( combi )
	dups1 <- combi[ dups ]
	dups_combi <- any( dups )
    PSF <- FALSE	
	#***** begin duplications of identifiers
    if ( dups_combi ){		
							
		NC <- max( table( table( combi) ) )		
		N1 <- nchar( paste(max(table(combi))))
		facets$psf <-  paste0("PF", 10^N1 + 1 )	
		for (cc in dups1){
			# cc <- dups1[1]
			ind <- which( combi == cc )
			N0 <- length(ind)
			h1 <- paste0("PF", 10^N1 + seq(1,N0) )	
			facets$psf[ind] <- h1
						}

		
		nullfacets <- c( nullfacets , "psf" )
        PSF <- TRUE
        cat("   -- Created a pseudo facet 'psf' (with zero effects)\n")
		cat("   -- because of non-unique person-facet combinations.\n") 		
					}
	#**** end duplications of identifiers
				
	# new formula
	formula_update <- paste( c( attr( tA , "term.labels") , nullfacets ) , collapse=" + ")
	formula_update <- stats::as.formula( paste0( "~ " , formula_update ) )
	xsi.setnull <- unique( c( xsi.setnull , nullfacets ) )
	
	if ( length(xsi.setnull)==0 ){
			xsi.setnull <- NULL
						}						

	#********************
	# dimensions for beta fixed
	D <- 1
	if ( ! is.null(B) ){
		D <- dim(B)[[3]]
					}
	if ( ! is.null(Q) ){
		D <- dim(Q)[[2]]
					}
	if ( is.null(beta.fixed) ){				
		beta.fixed <- cbind( 1 , 1:D , 0 )
							}

	res <- list( "formula_update" = formula_update , 
				"xsi.setnull" = xsi.setnull ,
				"beta.fixed" = beta.fixed ,
				"facets"=facets , "PSF" = PSF , "pid" = pid  )
	return(res)
		}
############################################################		