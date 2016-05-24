
#############################################################
# BIFIE.twolevelreg
BIFIE.twolevelreg <- function( BIFIEobj , dep , formula.fixed , formula.random ,
           idcluster , wgtlevel2 = NULL , wgtlevel1 = NULL , group=NULL , group_values=NULL , se=TRUE ,
		   globconv = 1E-6 , maxiter = 1000
				                  ){
				
				
		#****
		s1 <- Sys.time()
		cl <- match.call()			
		bifieobj <- BIFIEobj				

		#********	
		# extract variables						
		if (bifieobj$cdata){
			# formula_vars <- NULL
			# if (! is.null(formula) ){
			formula_vars <- c( base::all.vars( formula.fixed ) , 
									base::all.vars( formula.random ) )
			#					}
			varnames <- unique( c( dep ,  group , "one" , idcluster , formula_vars ,
								wgtlevel1 , wgtlevel2 ) )
			bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames=varnames )	
						}					
		# extract values		
		FF <- Nimp <- bifieobj$Nimp
		N <- bifieobj$N
		dat1 <- bifieobj$dat1
		wgt <- bifieobj$wgt
		wgtrep <- bifieobj$wgtrep
		varnames <- bifieobj$varnames
		RR <- bifieobj$RR
		datalistM <- bifieobj$datalistM
		fayfac <- bifieobj$fayfac	

		# create datalist
		datalistM <- as.data.frame( datalistM )
		colnames(datalistM) <- bifieobj$varnames
		bifieobj2 <- datalistM

		#*********** X predictor matrix
		m1 <- stats::model.matrix( formula.fixed , datalistM )
		m0 <- m1
		xnames <- colnames(m1)	

		m1 <- matrix( NA , nrow=nrow(bifieobj2) , ncol=ncol(m0) )
		m1[ match( rownames(m0),rownames(bifieobj2) ) , ] <- m0
		colnames(m1) <- colnames(m0)
		X_list <- as.matrix( m1 )

		#************* Z predictor matrix
		m1 <- stats::model.matrix( formula.random , datalistM )
		m0 <- m1
		znames <- colnames(m1)
		m1 <- matrix( NA , nrow=nrow(bifieobj2) , ncol=ncol(m0) )
		m1[ match( rownames(m0),rownames(bifieobj2) ) , ] <- m0
		colnames(m1) <- colnames(m0)
		Z_list <- as.matrix( m1 )

		#*************** y outcome values
		y_list <- as.vector( bifieobj2[ , dep ] )		
		globconv <- stats::var(y_list , na.rm=TRUE) * globconv		
		
		#*************** cluster identifiers
		
		dat1 <- bifieobj$dat1
		idcluster0 <- idcluster
		clusters <- unique( dat1[ , idcluster0 ] )
		idcluster <- match( dat1[ , idcluster0 ] , clusters )
		
		#************ weights		
		# weights
		wgttot <- wgt
		wgtlev2_full <- dat1[ , wgtlevel2 ]
		wgtlev2 <- stats::aggregate( wgtlev2_full , list(idcluster) , mean )[,2]
		if ( is.null(wgtlevel1) ){
			wgtlev1 <- wgttot / wgtlev2_full
			         } else {
			wgtlev1 <- dat1[ , wgtlevel1 ]
							}

		wgtrep1 <- wgtrep					
	
		#****** groups
		if ( is.null( group) ){ 
					nogroup <- TRUE 
								} else { 
					nogroup <- FALSE 
							}
		if (nogroup){
			group <- "one"
			group_values <- c(1)
				}		
		# group vector
		group_vec <- dat1[ , group ]

				
		if ( is.null(group_values ) ){ 		
			t1 <- fasttable( datalistM[ , group ] )				  
			group_values <- sort( as.numeric( paste( names(t1) ) ))
				}
		# use groups
		group_vec <- match( group_vec , group_values )
		group_values0 <- group_values
		group_values <- seq( 1 , length( unique( group_values) ) )
		
		group_vec <- as.numeric( group_vec)
		group_values <- as.numeric( group_values)
		
		
		
		if ( ! se ){ 
			wgtrep1 <- matrix( wgt , ncol=1 )
			RR <- 0
				}	

		
		#*********** estimate multilevel model		
		res <- .Call("bifie_mla2" , X_list , Z_list , y_list , wgttot , wgtlev2 , wgtlev1 ,
				 globconv , maxiter , group_vec - 1  , group_values - 1 , idcluster - 1 ,
				 wgtrep1 , Nimp , fayfac ,
				 PACKAGE="BIFIEsurvey")

				 
		# dimensions
		NX <- ncol(X_list)
		NZ <- ncol(Z_list)
		NC <- length(clusters)				
		GG <- length(group_values)	
		NP <- res$NP
		parnames <- paste0( "beta_" , xnames )
		for (ii in 1:NZ){
		  for (jj in 1:NZ){
		        if (ii==jj){
					v1 <- paste0( "Var_" , 	znames[ii] )
							}
		        if (ii<jj){
					v1 <- paste0( "Cov_" , 	znames[ii] , znames[jj] )
							}
		        if (ii>jj){
					v1 <- paste0( "Cor_" , 	znames[ii] , znames[jj] )
							}							
				parnames <- c( parnames , v1 )
							}
					}
		parnames <- c( parnames , "ResidVar")
		parnames <- c( parnames , "ExplVar_Lev2_Fixed" )
		parnames <- c( parnames , "ExplVar_Lev2_Random" )
		parnames <- c( parnames , "ResidVar_Lev2" )
		parnames <- c( parnames , "ExplVar_Lev1_Fixed" )
		parnames <- c( parnames , "ExplVar_Lev1_Random" )
		parnames <- c( parnames , "ResidVar_Lev1" )
		parnames <- c( parnames , "Var_Total" )
		parnames <- c( parnames , "R2_Lev2" )
		parnames <- c( parnames , "R2_Lev1" )
		parnames <- c( parnames , "R2_Total" )
		parnames <- c( parnames , "ICC_Uncond" )
		parnames <- c( parnames , "ICC_UncondWB" )
		parnames <- c( parnames , "ICC_Cond" )
				
		
		#*****
		# create parameter names
		
		
		#*******************
		# parameter table
		p1 <- parnames
		dfr <- data.frame( "parameter" = rep(p1,GG))
		ZZ <- NP
		# dfr$var <- rep(p2,GG)
		if (! nogroup){
		   dfr$groupvar <- group
		   dfr$groupval <- rep( group_values0 , each=ZZ )
					 }							 
		# dfr$Ncases <- rep( rowMeans( res$ncasesM ) , each=ZZ )
		# dfr$Nweight <- rep( rowMeans( res$sumwgtM ) , each=ZZ )	
					
#		dfr$est <- res$parsL$pars
#		dfr$SE <- res$parsL$pars_se
#		dfr$t <- round( dfr$est / dfr$SE , 2 )
#		dfr$p <- pnorm( - abs( dfr$t ) ) * 2
#		dfr$fmi <- res$parsL$pars_fmi
#		dfr$VarMI <- res$parsL$pars_varBetween
#		dfr$VarRep <- res$parsL$pars_varWithin

		dfr <- create_summary_table( res_pars=res$parsL , 
						 parsM=res$parsM   , parsrepM=res$parsrepM , 
						 dfr=dfr , BIFIEobj=BIFIEobj )				
		dfr <- clean_summary_table( dfr=dfr , RR=RR , se=se , Nimp=Nimp )
		
		
		
		# create vector of parameter names
		nogroupL <- rep( nogroup , nrow(dfr) )
		parnames <- paste0( dfr$parameter   , 
				ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
				ifelse( ! nogroupL , dfr$groupval , "" ) )
		rownames(dfr) <- parnames		
		
		#***************************
		# inference with mitools
		fvcovM <- res$fvcov
		vcov0 <- matrix(0 , nrow=GG*NX , ncol=GG*NX)
		vcov.list <- as.list(1:Nimp)
		for (ii in 1:Nimp){
			# ii <- 1
			vcov1 <- vcov0
			for (gg in 1:GG){
				# gg <- 1
				ind1 <- 1:NX + (gg-1)*NX
				vcov1[ ind1 , ind1 ] <- fvcovM[ ind1 , 1:NX + (ii-1)*NX ]
							}
			vcov.list[[ii]] <- vcov1
					}    

		parsM <- res$parsM
		NP <- nrow(parsM) / GG
		parsM0 <- matrix( 0 , nrow=NX*GG , ncol=Nimp)
		parnames_sel <- NULL
		for (gg in 1:GG){
			parsM0[ 1:NX + (gg-1)*NX , ] <- parsM[ 1:NX + (gg-1)*NP , ]
			parnames_sel <- c( parnames_sel , parnames[ 1:NX + (gg-1)*NP  ] )
                }
		parsM <- parsM0
		pars.list <- as.list(1:Nimp)
		for (ii in 1:Nimp){
		  pars.list[[ii]] <- parsM[,ii]
		  names(pars.list[[ii]] ) <- parnames_sel
						}
		if (Nimp==1){
			pars.list <- list( pars.list[[1]] , pars.list[[1]] )
			vcov.list <- list( vcov.list[[1]] , vcov.list[[1]] )
					}
		micombs <- mitools::MIcombine( pars.list , vcov.list )
		if ( ! se ){
			dfr$SE <- dfr$fmi <- dfr$VarRep <- NA
			v1 <- diag( micombs$variance )
			dfr[ parnames_sel , "SE" ] <- sqrt( v1 )
			dfr$t <- round( dfr$est / dfr$SE , 2 )
			dfr$p <- stats::pnorm( - abs( dfr$t ) ) * 2			
			dfr[ parnames_sel , "fmi" ] <- micombs$missinfo
			dfr$VarMI <- dfr$fmi * dfr$SE^2
			dfr$VarRep <- (1-dfr$fmi) * dfr$SE^2
				}
		
		#**************
		# variance decompositions
        vardecomp <- list( "Sigma_W_yXM" = res$Sigma_W_yXM , "Sigma_B_yXM" = res$Sigma_B_yXM ,
                "Sigma_W_yZM" = res$Sigma_W_yZM , "Sigma_B_yZM" = res$Sigma_B_yZM ,
				"totmean_yXM" = res$totmean_yXM , "totmean_yZM" = res$totmean_yZM )
				
		
		#*************************** OUTPUT ***************************************
		s2 <- Sys.time()
		timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
		res1 <- list( "stat" = dfr ,
				"output" = res , 
				"timediff" = timediff ,
				"Npers" = res$Npers , "Nclusters" = res$Nclusters , 
				"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
				"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
				"GG"=GG , "micombs" = micombs , "se" = se , 
				"parnames" = parnames , "parnames_sel" = parnames_sel ,
				"vardecomp" = vardecomp ,
				"CALL"= cl)		
		class(res1) <- "BIFIE.twolevelreg"		
		return(res1)																		
				}
#############################################################				


####################################################################################
# summary for BIFIE.linreg function
summary.BIFIE.twolevelreg <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat( paste0( "Number of persons:") , object$Npers , "\n")	
	cat( paste0( "Number of clusters:") , object$Nclusters , "\n\n")	
	
	cat("Statistical Inference for Two-Level Linear Regression \n\n")	
	obji <- object$stat
	rownames(obji) <- NULL
	print.object.summary( obji , digits=digits )
			}