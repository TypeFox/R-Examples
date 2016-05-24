

#######################################################################
# Path model
BIFIE.pathmodel <- function( BIFIEobj , lavaan.model , reliability=NULL , 
			group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	cl <- match.call()			
	bifieobj <- BIFIEobj
	dat <- bifieobj$dat1
	lavmodel <- lavaan.model

	# lavaanify syntax
	lavpartable <- TAM::lavaanify.IRT( lavmodel , dat )$lavpartable
	vars <- colnames(dat)
	
	lav1 <- lavpartable[ lavpartable$op %in% c("=~" , "~") , ]
	lav.vars <- unique( c( lav1$lhs , lav1$rhs ) )
	NV <- length(lav.vars)
	# observed variables
	obs.vars <- intersect( vars , lav.vars )
	# latent variables
	lat.vars <- setdiff( lav.vars , obs.vars )
	# create model matrix for latent variables
	NLV <- length(lat.vars)
	L <- matrix( 0 , nrow= NLV , ncol=NV )
	rownames(L) <- lat.vars
	colnames(L) <- lav.vars
	for (vv in 1:NLV){
		# vv <- 1
		lav.vv <- lav1[ ( lav1$op == "=~" ) & ( lav1$lhs == lat.vars[vv] ) , ]
		L[ vv , lav.vv$rhs ] <- ifelse( is.na(lav.vv$ustart) , 1 , lav.vv$ustart )
					}
									
	L_row_index <- match( rownames(L) , lav.vars )
					
	# model matrix for regressions
	lav2 <- lav1[ lav1$op == "~" , ]
	dep <- unique( lav2$lhs )
	NR <- length(dep)
	R <- matrix( 0 , nrow= NR , ncol=NV )
	rownames(R) <- dep
	colnames(R) <- lav.vars
	for (vv in 1:NR){
		lav.vv <- lav1[ ( lav1$op == "~" ) & ( lav1$lhs == dep[vv] ) , "rhs"]
		R[ vv , lav.vv ] <- 1
					}
	R_row_index <- match( rownames(R) , lav.vars )
	# matrix B with regression coefficients
	B <- matrix(0 ,nrow=NV , ncol=NV)
	rownames(B) <- colnames(B) <- lav.vars
	for (vv in 1:NR){
		B[ rownames(R)[vv] , ] <- R[vv , ]
					}
	# error variances and covariances
	E <- matrix( 0 , nrow=NV , ncol=NV )
	rownames(E) <- colnames(E) <- lav.vars
	lav2 <- lavpartable[ lavpartable$op == "~~" , ]
	lav2 <- stats::na.omit( lav2 )
	NG <- nrow(lav2)
	if (NG > 0 ){
	for ( gg in 1:NG){
		E[ lav2$lhs[gg] , lav2$rhs[gg] ] <- lav2$ustart[gg]
					}
				}
			
	# matrix power
	ind <- 0
	Bpow <- B
	Bpowsum <- Bpow
	for (kk in 1:100){
		Bpow <- Bpow %*% B
		Bpowsum <- Bpowsum + Bpow
		ind <- ind + 1
		if ( sum( Bpow ) < 1E-10 ){ break ; }
					}
	maxpow <- ind

	ind.vars <- colnames(B)[ which( colSums( B ) > 0 ) ]
	M1 <- matrix( ind.vars , ncol=1 )

	for (oo in 1:maxpow){
		NM1 <- nrow(M1)
		M2 <- NULL
		for (mm in 1:NM1){
			#mm <- 1
			v1 <- rownames(B)[  B[ , colnames(B) == paste(M1[ mm , oo ]) ] > 0  ]
			if (oo>1 ){ v1 <- c( "" , v1 ) }
			M1mm <-  t(M1[ mm , 1:oo ])
			HV <- length(v1)
			H1 <- matrix( "" , nrow=HV , ncol=oo+1 )
			H1[ , oo+1 ] <- v1
			H1[ , 1:oo ] <- matrix( M1mm , nrow=HV , ncol=oo, byrow=TRUE )
			M2 <- rbind(M2 , H1 )
					}							 
		M1 <- M2        
	}
	M1_index <- matrix( 0 , nrow(M1) , ncol(M1) )
	for (mm in 1:ncol(M1) ){
	   M1_index[ , mm ] <- match( M1[,mm] , lav.vars )
						}  
								   
	# compute total effects
	tot_paths <- t( apply( M1 , 1 , FUN = function(ll){
					vv <- ll[ ll != "" ]
					N1 <- length(vv)
					c( paste0( vv[1] , "->" , vv[N1] ) , N1 )
						} ) )
	tot_paths <- as.data.frame( tot_paths )
	tot_paths2 <- tot_paths[ as.numeric(paste(tot_paths[,2])) > 2 , ]
	paths <- unique( paste(tot_paths2[,1 ] ) )
	tot_paths$pathindex <- match( paste(tot_paths[,1]) , paths ) + nrow(M1)
	if ( mean( is.na(tot_paths$pathindex ) ) < 1 ){
		NP0 <- max( tot_paths$pathindex , na.rm=TRUE)
						} else { 
		NP0 <- nrow(M1) 
					}	
	
	
	
	coeff_index <- cbind( tot_paths$pathindex , M1_index )


	NV <- ncol(Bpowsum)
	indices <- which( Bpowsum > 0 )               

	
	if (bifieobj$cdata){
		varnames <- unique( c( obs.vars , group , "one" ) )
		bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames=varnames )	
						}				
	
	FF <- Nimp <- bifieobj$Nimp
	N <- bifieobj$N
	dat1 <- bifieobj$dat1
	wgt <- bifieobj$wgt
	wgtrep <- bifieobj$wgtrep
	varnames <- bifieobj$varnames
	RR <- bifieobj$RR
	datalistM <- bifieobj$datalistM
    fayfac <- bifieobj$fayfac	

	# vars_index <- match( lav.vars , colnames(dat) )
	vars_index <- match( lav.vars , colnames(dat1) )

	
	# unreliability
	unreliability <- rep(0 , NV)
	names(unreliability) <- lav.vars
	unreliability[ match( names(reliability) , lav.vars) ] <- 1 - reliability	
	NL <- nrow(L)
	if (NL == 0 ){
	    L <- matrix( 0 , nrow=1 , ncol=NV)
		colnames(L) <- lav.vars
				}
	
	tot_paths <- data.frame( tot_paths , "pathindex2" = tot_paths$pathindex )
	tot_paths$pathindex2 <- match( tot_paths$pathindex2 , 
			stats::na.omit(unique( tot_paths$pathindex2 )) ) + NP0
	tot_paths$pathindex2[ as.numeric(paste(tot_paths[,2])) <= 2 ] <- NA
	if ( mean( is.na(tot_paths$pathindex2 ) ) < 1 ){
		NP0 <- max( tot_paths$pathindex2 , na.rm=TRUE)
						} else { 
		NP0 <- nrow(M1) 
					}
	coeff_index1 <- cbind( coeff_index[,1] , tot_paths$pathindex2 , coeff_index[,-1] )
		
	wgt_ <- matrix( wgt , ncol=1 )
	if ( is.null( group) ){ nogroup <- TRUE } else { nogroup <- FALSE }
	cat(paste0( "|" , paste0( rep("*" , FF) , collapse="") , "|\n" ))
	if (nogroup){
	    group <- "one"
	    group_values <- c(1)
			}


	#@@@@***
    group_index <- match( group , varnames )
	#@@@@***

    if ( is.null(group_values ) ){ 
		t1 <- fasttable( datalistM[ , group_index ] )				  
	    group_values <- sort( as.numeric( paste( names(t1) ) ))
				}
	
	#@@@@***
	res00 <- BIFIE_create_pseudogroup( datalistM , group , group_index , group_values )				
	res00$datalistM -> datalistM 
	res00$group_index -> group_index
	res00$GR -> GR 
	res00$group_values -> group_values
	res00$group -> group
	#@@@@***			

					
						
	if (RR==1){ RR <- 0 }
	if ( ! se ){ 
		wgtrep <- matrix( wgt , ncol=1 )
		RR <- 0
				}							
																				
    #**** calculate path model
 #if (FALSE){	
	res <- .Call( "bifie_pathmodel" ,  datalistM , wgt_ , wgtrep ,
						vars_index - 1,  fayfac, Nimp , group_index - 1,
						group_values   , L , L_row_index - 1  , NL  , E ,
						R , R_row_index - 1 , coeff_index1 , NP0 , unreliability ,
						PACKAGE="BIFIEsurvey" ) 	

#	res <- bifie_pathmodel( datalistM , wgt_ , wgtrep ,	vars_index - 1,  fayfac, Nimp , group_index - 1,
#						group_values   , L , L_row_index - 1  , NL  , E , R , 
#                        R_row_index - 1 , coeff_index1 , NP0 , unreliability ) 
						
	GG <- length(group_values)

	
	#*** create parameter labels
	p0 <- paste(M1[,1])
	for (vv in seq(2,ncol(M1) ) ){
				p0 <- ifelse( M1[,vv] != "" , paste0( p0 , "->" , M1[,vv] ) , p0 ) 
						}
	
	p1 <- ifelse( tot_paths[,2] == "2" , paste0(M1[,2], "~" , M1[,1] ) , p0 )	

	p2 <- unique( paste(tot_paths[ ! is.na( tot_paths$pathindex ) , 1 ]) )
	p2 <- gsub( "->" , "-+>" , p2 )
	p1 <- c( p1 , p2 )
	p2 <- unique( paste(tot_paths[ ! is.na( tot_paths$pathindex2 ) , 1 ]) )
	p2 <- gsub( "->" , "-~>" , p2 )
	p1 <- c( p1 , p2 )	
	p1 <- c( p1 , paste0(p1 , "_stand") )
	p1 <- c( p1 , paste0( rownames(R) , "_R2" ) )
	p1 <- c( p1 , paste0( rownames(R) , "_ResidVar" ) )
  
	dfr <- data.frame( "parameter" = rep(p1,GG)	)
	# dfr$var <- rep(p1,GG)
	dfr$type <- ""
	ZZ <- length(p1) 
	
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( group_values , each=ZZ )
	             }		
			 
	dfr$Ncases <- rep( res$ncases[,1] , each=ZZ )
	dfr$Nweight <- rep( res$sumwgt[,1] , each=ZZ )	
	
	dfr <- create_summary_table( res_pars=res$parsL , 
				     parsM=res$parsM   , parsrepM=res$parsrepM , 
					 dfr=dfr , BIFIEobj=BIFIEobj )				
	dfr <- clean_summary_table( dfr=dfr , RR=RR , se=se , Nimp=Nimp )	
	
	dfr[ grep( "_R2" , paste(dfr$parameter) ) , "type"] <- "RSquared"
	dfr[ grep( "_ResidVar" , paste(dfr$parameter) ) , "type"] <- "ResidVar"		
	dfr[ grep( "~" , paste(dfr$parameter) ) , "type"] <- "RegrCoeff"			
	
	ind <- grep( "->" , paste(dfr$parameter) )
	if ( length(ind) > 0 ){
		dfr[ ind  , "type"] <- "PathCoeff"			
					}

	ind <- grep( "-+>" , paste(dfr$parameter) , fixed=TRUE )
	if ( length(ind) > 0 ){
		dfr[ ind  , "type"] <- "TotalEff"			
					}
	ind <- grep( "-~>" , paste(dfr$parameter) , fixed=TRUE )
	if ( length(ind) > 0 ){
		dfr[ ind  , "type"] <- "IndEff"			
					}					
					
	
	# create vector of parameter names
	nogroupL <- rep( nogroup , nrow(dfr) )
	parnames <- paste0( dfr$parameter   , "_" , 
			ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
			ifelse( ! nogroupL , dfr$groupval , "" ) )
	

	#@@@@***
	# multiple groupings
	dfr <- BIFIE_table_multiple_groupings( dfr , res00 )
	#@@@@***
					
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr ,
			"output" = res , 
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
			"GG"=GG , "parnames" = parnames , "CALL"= cl)
	class(res1) <- "BIFIE.pathmodel"
	return(res1)	
		}
###################################################################################

####################################################################################
# summary for BIFIE.pathmodel function
summary.BIFIE.pathmodel <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Statistical Inference for Path Model \n\n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}