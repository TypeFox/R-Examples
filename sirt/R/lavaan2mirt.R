
###################################################################
# Converting lavaan syntax into mirt syntax
lavaan2mirt <- function( dat , lavmodel , est.mirt =TRUE , poly.itemtype="gpcm" , ... ){
	# lavaanify model
	# lavmodel2 <- lavaan::lavaanify( lavmodel )		
#    lavmodel2 <- lavaanify.sirt( lavmodel )$lavpartable
	lavmodel2 <- TAM::lavaanify.IRT( lavmodel , data=dat )$lavpartable	
		
	# select used items
	items <- intersect( unique( paste(lavmodel2$rhs )) , colnames(dat) )
	dat <- dat[, items ]
	# maximum category
	maxK <- max( dat , na.rm=TRUE )
	ind1 <- which( lavmodel2$op == "~~")
	ind2 <- which( ( lavmodel2$lhs == lavmodel2$rhs ) & ( lavmodel2$lhs %in% items ) )
	ind <- intersect(ind1 , ind2 )
	if ( length(ind) > 0 ){
		lavmodel2 <- lavmodel2[ - ind , ]
						}
						
	# variable names
	items <- colnames(dat)
	# extract factors
	factors <- unique( paste(lavmodel2$lhs[ lavmodel2$op == "=~" ] ))
	# create mirt syntax
	mirtmodel <- NULL
	sel.items <- NULL
	#------------
	#**** loop over factors
	for (ff in factors){
		# ff <- factors[1]
		ind.ff <- which( ( lavmodel2$op == "=~"   ) & ( lavmodel2$lhs == ff ) )
		mff <- match( lavmodel2[ ind.ff , "rhs" ] , items  )
		sel.items <- c( sel.items , mff )
		mirtmodel <- paste0( mirtmodel , 
				paste0( ff , "=" , paste0( mff , collapse="," ) , "\n")  )
						}				
	#------------						
	#***** look for constraints	among item parameters (loadings and intercepts)				
	lavmodel21 <- lavmodel2
	lavmodel21 <- lavmodel21[ paste(lavmodel21$op) %in% c("=~","|","?=") , ]	
	lavlabels <- unique(paste(lavmodel21$label))	
	lavlabels <- lavlabels[ paste(lavlabels ) != "" ]
	
	if ( length(lavlabels) > 0 ){
		vv0 <- "CONSTRAIN = "
		for (ll in lavlabels ){
			# ll <- lavlabels[1]
			lav2.ll <- lavmodel21[ paste(lavmodel21$label) == ll , ]
			if (lav2.ll$op[1] == "=~"){ 
				pars.ll <- paste0( "a" , match( lav2.ll$lhs , factors )[1] )
				isel <- lav2.ll$rhs
								}
			# Example: CONSTRAIN = (1-12,a1),(1-12,a2),(1-12,a3)
			if (lav2.ll$op[1] == "|"){ 
				pars.ll <- gsub( "t" , "d" , lav2.ll$rhs[1] )
				if (maxK==1){ pars.ll <- "d" }
				# if (maxK>1){ pars.ll <- gsub( "|" , "" , pars.ll ) }				
				isel <- lav2.ll$lhs
								}
			if (lav2.ll$op[1] == "?="){ 
				pars.ll <- gsub( "g1" , "g" , lav2.ll$rhs[1] )
				pars.ll <- gsub( "s1" , "u" , pars.ll )				
				# if (maxK>1){ pars.ll <- gsub( "|" , "" , pars.ll ) }				
				isel <- lav2.ll$lhs
								}


								
			vv2 <- paste0( paste0( match( isel , items ) , collapse=",") , "," , pars.ll )
			vv2 <- paste0("(" , vv2 , ")")
			if (ll != lavlabels[1] ){ vv0 <- paste0( vv0 , "," ) }
			vv0 <- paste0( vv0 , vv2 )
							}
		mirtmodel <- paste0( mirtmodel , vv0  , "\n")		
						}

	#------------
	#**** estimate variances and covariances
	lavmodel21 <- lavmodel2
	lavmodel21 <- lavmodel21[ paste(lavmodel21$op) %in% c("~~") , ]							
	LL <- nrow(lavmodel21)
	vv1 <- "COV = "
	if (LL>0){
		for (ll in 1:LL){
			 vv2 <- paste0(lavmodel21[ll,"lhs"] , "*" , lavmodel21[ll,"rhs"] )
			 if (ll>1){ vv1 <- paste0( vv1 , "," ) }
			 vv1 <- paste0( vv1 , vv2 )	
					}
		mirtmodel <- paste0( mirtmodel , vv1 , "\n")
				}
	#------------
	#**** estimate means
	lavmodel21 <- lavmodel2
	lavmodel21 <- lavmodel21[ paste(lavmodel21$op) %in% c("~1") , ]							
	lavmodel21 <- lavmodel21[ paste(lavmodel21$lhs) %in% factors , ]	
	LL <- nrow(lavmodel21)
	vv1 <- "MEAN = "
	if (LL>0){
		for (ll in 1:LL){
		     vv2 <- paste0(lavmodel21[ll,"lhs"])
			 if (ll>1){ vv1 <- paste0( vv1 , "," ) }
			 vv1 <- paste0( vv1 , vv2 )	
					}
		mirtmodel <- paste0( mirtmodel , vv1 , "\n")
				}				
	#------------						
	#**** create object of class mirt.model
	mirtmodel1 <- mirt::mirt.model(mirtmodel)
	#------------	
	#**** create parameter values

	# define item types: 4PL or gpcm
	I <- ncol(dat)
	maxK1 <- apply( dat , 2 , max , na.rm=TRUE )

	typeK1 <- "4PL"
	ind1 <- which( lavmodel2$op == "=~" )
	ind2 <- which( ( ! lavmodel2$free  ) & ( lavmodel2$ustart %in% c(0,1) ) )
	ind2 <- intersect( ind1 , ind2 )
    if ( length(ind1) == length(ind2) ){ typeK1 <- "Rasch" }	

	itemtype <- ifelse( maxK1==1 , typeK1 ,  poly.itemtype )	
	
	mirtpars <- mirt::mirt( dat , mirtmodel1 , itemtype=itemtype , pars="values")
		
	#***** handle guessing parameters
	lavmodel21 <- lavmodel2[ lavmodel2$op == "?=" , ]
	lavmodel21 <- lavmodel21[ lavmodel21$rhs == "g1" , ]
	lavmodel21 <- lavmodel21[ lavmodel21$free > 0 , ]

	
	if ( nrow(lavmodel21) > 0 ){	
			ind <- which( ( mirtpars$item %in% lavmodel21$lhs ) &
							( mirtpars$name == "g") )
			ind <- match( paste0(lavmodel21$lhs,"-","g") , 
						paste0(mirtpars$item,"-",mirtpars$name) )													
			mirtpars[ind,"value"] <- .25
			mirtpars[ind,"est"] <- TRUE			
						} else {	# no guessing parameters
			ind <- which( ( mirtpars$item %in% items ) &
							( mirtpars$name == "g") )
			if ( length( ind ) > 0 ){			
				mirtpars[ind,"value"] <- 0
				mirtpars[ind,"est"] <- FALSE			
									}
							}
	#***** handle slipping parameters
	lavmodel21 <- lavmodel2[ lavmodel2$op == "?=" , ]
	lavmodel21 <- lavmodel21[ lavmodel21$rhs == "s1" , ]
	lavmodel21 <- lavmodel21[ lavmodel21$free > 0 , ]
	
	
	if ( nrow(lavmodel21) > 0 ){	
			ind <- match( paste0(lavmodel21$lhs,"-","u") , 
						paste0(mirtpars$item,"-",mirtpars$name) )
			mirtpars[ind,"value"] <- .95
			mirtpars[ind,"est"] <- TRUE			
						}  else {	# no guessing parameters
			ind <- which( ( mirtpars$item %in% items ) &
							( mirtpars$name == "u") )
			if ( length( ind ) > 0 ){							
				mirtpars[ind,"value"] <- 1
				mirtpars[ind,"est"] <- FALSE			
								}
							}
	
	#------------	
	#**** constraints for parameter values
    lavmodel21 <- lavmodel2[ lavmodel2$free == 0 , ]
	LL <- nrow(lavmodel21)
	if (LL>0){
		for (ll in 1:LL){
	#		ll <- 1
	        item.ll <- NULL
			lavmodel21.ll <- lavmodel21[ll,]
			# factor loadings
			if ( lavmodel21.ll$op == "=~" ){
				par.ll <- paste0("a",match( paste0(lavmodel21.ll$lhs) , factors ) )
				item.ll <- paste0(lavmodel21.ll$rhs)
										}
			# thresholds
			if ( lavmodel21.ll$op == "|" ){
				par.ll <- gsub( "t" , "d" , lavmodel21.ll$rhs[1] )
				if (maxK==1){ par.ll <- "d" }
				item.ll <- paste0(lavmodel21.ll$lhs)
										}
			# guessing/slipping parameters
			if ( lavmodel21.ll$op == "?=" ){
				par.ll <- gsub( "g1" , "g" , lavmodel21.ll$rhs[1] )
				par.ll <- gsub( "s1" , "u" , par.ll )				
				if (par.ll=="u"){ 
						lavmodel21.ll$ustart <- 1 - lavmodel21.ll$ustart  
								}			
				item.ll <- paste0(lavmodel21.ll$lhs)
										}									
																				
			# covariances
			if ( lavmodel21.ll$lhs %in% colnames(dat) ){
					lavmodel21.ll$op <- "" 
						}
			if ( lavmodel21.ll$op == "~~" ){
			    item.ll <- "GROUP"
				par.ll <- c( match( lavmodel21.ll$lhs , factors ) ,
					match( lavmodel21.ll$rhs , factors ) )
				par.ll <- paste0("COV_" ,
						paste0( sort(par.ll, decreasing=TRUE ) , collapse="" ) )
										}
			# covariances
			if ( lavmodel21.ll$op == "~1" ){
			    item.ll <- "GROUP"
				par.ll <- match( lavmodel21.ll$lhs , factors ) 
				if ( is.na( par.ll ) ){ item.ll <- NULL }				
				par.ll <- paste0("MEAN_" , par.ll )			
										}			
			#++ parameter fixing
			if ( ! is.null(item.ll) ){										
				ind.ll <- which( ( mirtpars$item == item.ll ) & ( mirtpars$name == par.ll ) )	
				mirtpars[ind.ll,"est"] <- FALSE
				mirtpars[ind.ll,"value"] <- lavmodel21.ll$ustart
								}
							}
				}

	#------------
	#****	
    # estimate mirt model
	res <- NULL

	if ( est.mirt ){
		res$mirt <- mirt::mirt( dat , model=mirtmodel1 , pars=mirtpars , 
						itemtype=itemtype ,  ... )
					}   # else {
		res$mirt.model <- mirtmodel1				
		res$mirt.syntax <- mirtmodel
		res$mirt.pars <- mirtpars
		res$lavaan.model <- lavmodel2
		res$dat <- dat  
#				}
    return(res)
	}
###################################################################