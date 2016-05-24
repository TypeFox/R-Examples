

#######################################################################
# cross tabulation
BIFIE.crosstab <- function( BIFIEobj , vars1 , vars2 , 
		    vars_values1 = NULL , vars_values2 = NULL , 
			group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	cl <- match.call()	
	bifieobj <- BIFIEobj
	vars1 <- vars1[1]
	vars2 <- vars2[1]	
	if (bifieobj$cdata){
		varnames <- unique( c( vars1 , vars2 , group , "one") )
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

	
	if (RR==1){ RR <- 0 }
	if ( ! se ){ 
		wgtrep <- matrix( wgt , ncol=1 )
		RR <- 0
				}	
	

	vars_index1 <- which( varnames == vars1 )
	vars_index2 <- which( varnames == vars2 )
    # vars values
	if ( is.null(vars_values1 ) ){ 
		# t1 <- table( dat1[ , vars_index1 ] )	
		t1 <- fasttable( datalistM[ , vars_index1 ] )		
	    vars_values1 <- sort( as.numeric( paste( names(t1) ) ))
				}
	if ( is.null(vars_values2 ) ){ 
		# t1 <- table( dat1[ , vars_index2 ] )				  
		t1 <- fasttable( datalistM[ , vars_index2 ] )		
	    vars_values2 <- sort( as.numeric( paste( names(t1) ) ))
				}
	
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
				
				
	GG <- length(group_values)
				
					
	#**************************************************************************#
	# Rcpp call
	res <- .Call( "bifie_crosstab" , datalistM , wgt_ , wgtrep , vars_values1  ,
		vars_index1 - 1,    vars_values2 ,  vars_index2 - 1 , fayfac ,
		Nimp , group_index - 1 , group_values , PACKAGE="BIFIEsurvey")
	ZZ <- nrow(res$ctparsM)
	design_pars <- res$design_pars
	VV1 <- length(vars_values1)
	VV2 <- length(vars_values2)
	#*********
	# joint distributions
	dfr1 <- data.frame("var1" = vars1[1] , "varval1" = design_pars[,1] )
	dfr1$var2 <- vars2[1]
	dfr1$varval2 <- design_pars[,2]
	dfr1$group <- group
	dfr1$groupval <- design_pars[,3]
	XX1 <- nrow(dfr1)
	dfr1$Ncases <- rowMeans( res$ncasesM )
	dfr1$Nweight <- rowMeans( res$sumwgtM )	
	XX2 <- 3*XX1
	## // probs_joint    ZZ
	## // probs_rowcond  ZZ
	## // probs_colcond   ZZ
	## // probs_rowmarg   VV1*GG
	## // probs_colmarg  VV2*GG
	dfr1 <- data.frame("prob" =  rep( c("joint" , "rowcond" , "colcond") , each=XX1 ) ,
	          dfr1[ rep(1:XX1 , 3 ) , ] )
	dfr1$est <- res$ctparsL$pars[ 1:XX2 ]
	dfr1$SE <- res$ctparsL$pars_se[ 1:XX2 ]		
	dfr1$fmi <- res$ctparsL$pars_fmi[ 1:XX2 ]
	dfr1$df <- rubin_calc_df( res$ctparsL , Nimp , indices = 1:XX2 )	
	dfr1$VarMI <- res$ctparsL$pars_varBetween[ 1:XX2 ]
	dfr1$VarRep <- res$ctparsL$pars_varWithin[ 1:XX2 ]
	rownames(dfr1) <- NULL		  
	parnames <- paste0( dfr1$prob , "_" , dfr1$var1 , dfr1$val1 , "_" ,
				dfr1$var2 , dfr1$val2 , "_" , dfr1$group , dfr1$groupval )
		if (BIFIEobj$NMI ){
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$ctparsM[1:XX2,] , parsrepM=res$ctparsrepM[1:XX2,] , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )			
			dfr1$est <- res1$pars
			dfr1$SE <- res1$pars_se
			# dfr$t <- round( dfr$perc / dfr$perc_SE , 2 )
			dfr1$df <- res1$df
			# dfr$p <- pt( - abs( dfr$t ) , df=dfr$df) * 2			
			dfr1$fmi <- res1$pars_fmi
		    dfr1$fmi_St1 <- res1$pars_fmiB
		    dfr1$fmi_St2 <- res1$pars_fmiW						
			dfr1$VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr1$VarMI_St1 <- res1$pars_varBetween1
			dfr1$VarMI_St2 <- res1$pars_varBetween2						
			dfr1$VarRep <- res1$pars_varWithin	
							}
				
	#*****
	# marginal distributions
	XX3 <- GG*(VV1+VV2)	
	dfr2 <- data.frame( "prob" = c( rep( "rowmarg" , VV1*GG) , rep( "colmarg" , VV2*GG) ) )
	dfr2$var <- c( rep(vars1,VV1*GG) , rep(vars2,VV2*GG) )	
	dfr2$varval <- c( rep( vars_values1 , GG ) , rep( vars_values2 , GG ) )
	dfr2$group <- group
	dfr2$groupval <- c( rep( group_values , each = VV1 ) , rep( group_values , each = VV2 )  )
	l1 <- seq( XX2+1 , XX2 + XX3 )
	dfr2$est <- res$ctparsL$pars[ l1 ]
	dfr2$SE <- res$ctparsL$pars_se[ l1 ]
	dfr2$fmi <- res$ctparsL$pars_fmi[ l1 ]
	dfr2$df <- rubin_calc_df( res$ctparsL , Nimp , indices = l1)
	dfr2$VarMI <- res$ctparsL$pars_varBetween[ l1 ]
	dfr2$VarRep <- res$ctparsL$pars_varWithin[ l1 ]
	parnames2 <- paste0( dfr2$prob , "_" , dfr2$var , dfr2$val , "_" ,
				 dfr2$group , dfr2$groupval )
	parnames <- c( parnames , parnames2 )

		if (BIFIEobj$NMI ){
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$ctparsM[ l1 ,] , parsrepM=res$ctparsrepM[ l1 ,] , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )			
			dfr2$est <- res1$pars
			dfr2$SE <- res1$pars_se
			# dfr$t <- round( dfr$perc / dfr$perc_SE , 2 )
			dfr2$df <- res1$df
			# dfr$p <- stats::pt( - abs( dfr$t ) , df=dfr$df) * 2			
			dfr2$fmi <- res1$pars_fmi
		    dfr2$fmi_St1 <- res1$pars_fmiB
		    dfr2$fmi_St2 <- res1$pars_fmiW						
			dfr2$VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr2$VarMI_St1 <- res1$pars_varBetween1
			dfr2$VarMI_St2 <- res1$pars_varBetween2						
			dfr2$VarRep <- res1$pars_varWithin	
							}
	
	
	#*****
	# effect sizes
	## // w_es  2*GG
	## // gamma_es  GG
	## // lambda   3*GG
	## // kruskal_tau  3*GG	
	XX4 <- nrow(dfr1) + nrow(dfr2)	
	XX5 <- (2+1+3+3)*GG
	dfr3 <- data.frame( "parm" = c( rep("w",GG) , rep("V",GG) , 
	          rep("gamma",GG) ,
			rep(c("lambda", "lambda_X" ,"lambda_Y"),GG) , 
			rep( c("tau","tau_X","tau_Y"), GG ) ) )
	dfr3$group <- group
	dfr3$groupval <- c(rep(group_values ,1) , rep(group_values, 1) , rep(group_values,1) ,
				rep(group_values,each=3) , rep(group_values,each=3) )
	l1 <- seq( XX4+1 , XX4 + XX5 )
	dfr3$est <- res$ctparsL$pars[ l1 ]
	dfr3$SE <- res$ctparsL$pars_se[ l1 ]
	dfr3$fmi <- res$ctparsL$pars_fmi[ l1 ]
	dfr3$df <- rubin_calc_df( res$ctparsL , Nimp , indices = l1)
	dfr3$VarMI <- res$ctparsL$pars_varBetween[ l1 ]
	dfr3$VarRep <- res$ctparsL$pars_varWithin[ l1 ]	
	
		if (BIFIEobj$NMI ){
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$ctparsM[ l1 ,] , parsrepM=res$ctparsrepM[ l1 ,] , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )			
			dfr3$est <- res1$pars
			dfr3$SE <- res1$pars_se
			# dfr$t <- round( dfr$perc / dfr$perc_SE , 2 )
			dfr3$df <- res1$df
			# dfr$p <- stats::pt( - abs( dfr$t ) , df=dfr$df) * 2			
			dfr3$fmi <- res1$pars_fmi
		    dfr3$fmi_St1 <- res1$pars_fmiB
		    dfr3$fmi_St2 <- res1$pars_fmiW						
			dfr3$VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr3$VarMI_St1 <- res1$pars_varBetween1
			dfr3$VarMI_St2 <- res1$pars_varBetween2						
			dfr3$VarRep <- res1$pars_varWithin	
							}	
	
	
	parnames3 <- paste0( dfr3$parm , "_" ,
				 dfr3$group , dfr3$groupval )
	parnames <- c( parnames , parnames3 )
	if ( ( ! se ) &  ( RR==0 ) ){				
		dfr1$df <- dfr1$SE <- dfr1$fmi <- dfr1$VarMI <- dfr1$VarRep <- NULL
		dfr2$df <- dfr2$SE <- dfr2$fmi <- dfr2$VarMI <- dfr2$VarRep <- NULL		
		dfr3$df <- dfr3$SE <- dfr3$fmi <- dfr3$VarMI <- dfr3$VarRep <- NULL				
				}				

	if ( Nimp == 1 ){				
		dfr1$fmi <- dfr1$VarMI <-NULL
		dfr2$fmi <- dfr2$VarMI <-  NULL		
		dfr3$fmi <- dfr3$VarMI <-  NULL				
				}				
				
	# create vector of parameter names
#	nogroupL <- rep( nogroup , nrow(dfr) )
#	parnames <- paste0( dfr$var   , "_" , dfr$varval , 
#			ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
#			ifelse( ! nogroupL , dfr$groupval , "" ) )	
# parnames <- NULL	
	
	# compute ad hoc chi square statistics (without resampling)
	ncases_gg <- res$ncases_ggM	
	l1 <- seq( XX4+1 , XX4 + GG)
	wes <- res$ctparsM[l1,]
	chisquare <- wes^2 * ncases_gg
	p_chi2 <- (VV1-1)*(VV2-1)
	p_chi2 <- 1- stats::pchisq( chisquare , df=p_chi2 )
	dfr4 <- data.frame("group" = group , "groupval" = group_values )
    for (ii in 1:GG){
		m1 <- miceadds::micombine.chisquare( dk = chisquare[ii,] , df = (VV1-1)*(VV2-1) , 
					display=FALSE )
		dfr4[ii,"chi2"] <- m1["D"]
		dfr4[ii , "df" ] <- m1["df"] 
		dfr4[ ii , "p"] <- m1["p"]
				}		
				

	#@@@@***
	# multiple groupings
	dfr1 <- BIFIE_table_multiple_groupings( dfr1 , res00 )
	#@@@@***

	#@@@@***
	# multiple groupings
	dfr2 <- BIFIE_table_multiple_groupings( dfr2 , res00 )
	#@@@@***
						
				
				
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat.probs" = dfr1 , "stat.marg" = dfr2 , 
			"stat.es" = dfr3 , "chisq.test" = dfr4 ,
			"output" = res , "timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
			"parnames" = parnames , "CALL"= cl )
	class(res1) <- "BIFIE.crosstab"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.crosstab function
summary.BIFIE.crosstab <- function( object , digits=3 , ... ){
    BIFIE.summary(object)
	cat("Joint and Conditional Probabilities\n")	
	obji <- object$stat.probs
	print.object.summary( obji , digits=digits )
	cat("\nMarginal Probabilities\n")	
	obji <- object$stat.marg
	print.object.summary( obji , digits=digits )	
	cat("\nEffect Sizes\n")	
	obji <- object$stat.es
	print.object.summary( obji , digits=digits )		
	cat("\nChi Squared Test (Ignoring Weights and Sampling Scheme)\n")	
	obji <- object$chisq.test
	print.object.summary( obji , digits=digits )		
			}