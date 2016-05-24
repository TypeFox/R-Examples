
#################################################################
# transforming data in BIFIE.data object
BIFIE.data.transform <- function( bifieobj , transform.formula ,
     varnames.new = NULL ){
	varnames <- bifieobj$varnames	
	transform.formula <- stats::as.formula( transform.formula )
	# select variables which should be transformed
	vars <- base::all.vars( transform.formula )
	# ind_vars <- which( varnames %in% vars )	
	ind_vars <- unlist( sapply( vars , FUN = function(vv){
	             which( varnames == vv )
							} ) )		
	# check whether all variables are included in the BIFIEdata object
	h1 <- setdiff( vars , varnames )
	if ( length(h1)>0 ){
	    stop( "Following variables are not in the BIFIEdata object: \n\n   " , 
				paste( h1 , collapse=" " ) )
					}		
	#***--- distinction BIFIEdata and BIFIEcdata
	cdata <- bifieobj$cdata
	if ( ! cdata ){
		dfr <- base::as.data.frame( bifieobj$datalistM[ , ind_vars] )
		colnames(dfr) <- vars
				 }
	if ( cdata ){
		dfr1 <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames = vars )
		dfr <- as.data.frame( dfr1$datalistM )
		colnames(dfr) <- dfr1$varnames
				 }
	N <- bifieobj$N
	N1 <- bifieobj$Nimp * N
	N2 <- ncol( bifieobj$dat1)
	dfr_long <- dfr			 
	Nimp <- bifieobj$Nimp
	#****
	# check whether some variables should be removed in original BIFIE.data object
	if ( ! is.null( varnames.new) ){
		varnames.old <- bifieobj$varnames
		select_vars <- setdiff( varnames.old , varnames.new )
        bifieobj <- BIFIEdata.select( 	bifieobj , varnames = select_vars )	
		
		# removed variables
		rm_vars <- intersect( varnames.old , varnames.new )
		if ( length(rm_vars) > 0 ){
			cat( paste0("Removed "  , length(rm_vars) , 
						" original variables: " , paste0( rm_vars , collapse = " " ) , "\n") )
								 }
		varnames <- bifieobj$varnames								 
						}


	#***---	
	#*** construction of new variables
	M1_long <- NULL
	for (ii in 1:Nimp){
		# ii <- 1
		dfr <- as.data.frame( dfr_long[ ( ii-1)*N + 1:N , ] )
		colnames(dfr) <- colnames(dfr_long)
		rownames(dfr) <- 1:N		
		M1 <- model.matrix( transform.formula , dfr )
		varnames.added <- colnames(M1)
		varsindex.added <- seq( N2 + 1 , N2 + ncol(M1) )
		M1.new <- matrix( NA , nrow=N , ncol=ncol(M1) )
		colnames(M1.new) <- varnames.added
		VV <- ncol(M1)
		M1.new[ match( rownames(M1),rownames(dfr) ) , ] <- M1			
		M1_long <- rbind( M1_long , M1.new )
					}
	# resulting object is M1.new
	M1.new <- M1_long	
	
	#***---		
	varnames.added1 <- varnames.added
	if ( ! is.null(varnames.new) ){
		V21 <- length(varnames.added)
        V22 <- length(varnames.new)
        varnames.added[ seq(1 , min(V21,V22)) ] <- varnames.new[ seq(1,min(V21,V22) ) ]
				}
	varnames1 <- c( varnames , varnames.added )		
	
	#***--- distinction BIFIEdata and BIFIEcdata
	if ( ! cdata ){
		bifieobj$datalistM <- as.matrix( cbind( bifieobj$datalistM , M1.new ) )
		colnames(bifieobj$datalistM) <- NULL		
    	bifieobj$dat1 <- as.matrix( bifieobj$datalistM[ seq( N*(Nimp-1) + 1 , Nimp*N ) , ,drop=FALSE])				
		colnames(bifieobj$dat1) <- varnames1		
				}
						
	if ( cdata ){			
		M1.new <- as.matrix(M1.new)				
		VV2 <- ncol(bifieobj$dat1)	
		# create indicators
        res2 <- .Call( "bifie_bifiedata2bifiecdata" ,
					       M1.new , bifieobj$Nimp ,  
						   PACKAGE="BIFIEsurvey" )												
					   
		# colnames(res2$datalistM_imputed) <- c("_imp" , "subj" , "variable" , "value")						   
		datalistM_ind <- res2$datalistM_ind
		datalistM_imputed <- res2$datalistM_imputed	
	    datalistM_impindex <- res2$datalistM_impindex		

#        datalistM_imputed[,"variable"] <- datalistM_imputed[,"variable"] + VV2							
        datalistM_impindex[,2] <- datalistM_impindex[,2] + VV2										
		bifieobj$datalistM_imputed <- rbind( bifieobj$datalistM_imputed , datalistM_imputed )
#		bifieobj$datalistM_ind <- cbind( bifieobj$datalistM_ind , datalistM_ind )
    	dat1 <- as.matrix( M1.new[ seq( N*(Nimp-1) + 1 , Nimp*N ) , ])
		bifieobj$dat1 <- cbind( bifieobj$dat1 , dat1 )		
		colnames(bifieobj$dat1) <- varnames1
		bifieobj$datalistM_impindex <- rbind( bifieobj$datalistM_impindex , datalistM_impindex )
		bifieobj$datalistM_ind <- cbind( bifieobj$datalistM_ind , datalistM_ind )
				}

						
	#*****
	# include variable names
	bifieobj$varnames <- varnames1
	bifieobj$varnames.added <- varnames.added
	bifieobj$varsindex.added <- varsindex.added
	cat( paste0( "Included " , VV , " variables: " , paste0( varnames.added , collapse=" " ) , "\n") )
	#**** add variable names in list
	dfr3 <- bifieobj$variables
	VV2 <- length(bifieobj$varnames.added)
    n0 <- max( dfr3$index ) 
	dfr3[  , "variable"]  <- bifieobj$varnames[ seq( 1, nrow(dfr3) ) ]
	dfr2 <- data.frame( "index" = n0 + 1:VV2 , 
				"variable" = varnames.added ,
				"variable_orig" = varnames.added1 , 
				"source"= paste0( as.character(transform.formula) , collapse=" ") )	
	dfr3 <- rbind( dfr3 , dfr2 )
	dfr3 -> bifieobj$variables
	bifieobj$Nvars <- ncol(bifieobj$dat1)
	return( bifieobj )
		}
#################################################################		