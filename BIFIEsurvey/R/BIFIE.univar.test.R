

#######################################################################
# BIFIE univar.test
BIFIE.univar.test <- function( BIFIE.method , wald_test=TRUE ){
	#****
	s1 <- Sys.time()	
	cl <- match.call()		
	res5 <- BIFIE.method	
	if (res5$group == "one"){ stop("This function can only be applied with a grouping variable.\n")}
	
	if (BIFIE.method$RR < 2){ wald_test <- FALSE }
	
	mean1M <- res5$output$mean1M
	sd1M <- res5$output$sd1M
	sumweightM <- res5$output$sumweightM
	GG <- res5$GG
	group_values <- ( res5$stat$groupval )[1:GG]
	if( res5$group == "pseudogroup" ){
		is_pseudogroup <- TRUE
						} else {
		is_pseudogroup <- FALSE				
						}						
	if ( is.null( group_values) ){
		group_values <- 1:GG
								}
	
	mean1repM <- res5$output$mean1repM
	sd1repM <- res5$output$sd1repM
	sumweightrepM <- res5$output$sumweightrepM
	fayfac <- res5$fayfac
	VV <- res5$VV
	vars <- res5$vars
    RR <- res5$RR
	if (RR==1){ RR <- 0 }
	group <- res5$group
	N <- res5$N
	Nimp <- res5$Nimp

	
	#*****
	# Rcpp call
	res <- .Call("bifie_test_univar" , mean1M , sd1M , sumweightM , GG , group_values ,
			mean1repM , sd1repM , sumweightrepM  , fayfac , PACKAGE="BIFIEsurvey")
#	res <- bifie_test_univar(  mean1M , sd1M , sumweightM , GG , group_values ,
#			mean1repM , sd1repM , sumweightrepM  , fayfac )

			
	#****
	# output eta^2		
	dfr <- data.frame( "var" = vars	, "group" = group  )	
    dfr$eta2 <- res$eta2L$pars^2	
	dfr$eta <- res$eta2L$pars
	dfr$eta_SE <- res$eta2L$pars_se
	dfr$fmi <- res$eta2L$pars_fmi
	dfr$df <- rubin_calc_df( res$eta2L , Nimp )		
	dfr$VarMI <- res$eta2L$pars_varBetween
	dfr$VarRep <- res$eta2L$pars_varWithin	
		if (BIFIE.method$NMI ){
			# eta2
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$eta2M , parsrepM=res$eta2repM , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIE.method$Nimp_NMI , comp_cov = FALSE )			
			dfr$eta <- res1$pars
			dfr$eta_SE <- res1$pars_se
			dfr$df <- res1$df
			dfr$eta_fmi <- res1$pars_fmi
			dfr$eta_VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr$eta_VarRep <- res1$pars_varWithin		
							}	
	if (RR==0){				
		dfr$df <- dfr$SE <- dfr$fmi <- dfr$VarMI <- dfr$VarRep <- NULL
				}				
	# extract replicated statistics for d and eta squared							
	stat.eta2 <- dfr
	
	#****
	# output d values
	group_values_matrix <- res$group_values_matrix	
	ZZ <- nrow(group_values_matrix)
	dfr <- data.frame( "var" = rep(vars,each=ZZ)	)	
	g2 <- group_values_matrix[ rep(1:ZZ , VV) , , drop=FALSE]
	g2 <- data.frame(g2)
	colnames(g2) <- c("groupval1" , "groupval2")
	dfr$group <- group
	dfr <- cbind( dfr , g2 )
	r5 <- res5$stat
	h1 <- NULL
	h4 <- h3 <- h2 <- NULL
    for (kk in 1:2){
		group_values_matrix[,kk] <- match( group_values_matrix[,kk] , group_values )
		}
	for (vv in 1:VV){
		h1 <- c(h1 , r5[ (vv-1)*GG + group_values_matrix[,1] , "M" ] )
		h2 <- c(h2 , r5[ (vv-1)*GG + group_values_matrix[,2] , "M" ] )		
		h3 <- c(h3, r5[ (vv-1)*GG + group_values_matrix[,1] , "SD" ] )
		h4 <- c(h4 , r5[ (vv-1)*GG + group_values_matrix[,2] , "SD" ] )		
					}
	dfr$M1 <- h1
    dfr$M2 <- h2
    dfr$SD <- sqrt( ( h3^2 + h4^2 ) / 2 )		
	dfr$d <- res$dstatL$pars
	dfr$d_SE <- res$dstatL$pars_se
	dfr$d_t <- dfr$d / dfr$d_SE
	dfr$d_df <- rubin_calc_df( res$dstatL , Nimp )	
	dfr$d_p <- stats::pt(  - abs(dfr$d_t ) , df = dfr$d_df)*2
	dfr$d_fmi <- res$dstatL$pars_fmi
	dfr$d_VarMI <- res$dstatL$pars_varBetween
	dfr$d_VarRep <- res$dstatL$pars_varWithin


	
	if ( BIFIE.method$NMI ){
		res1 <- BIFIE_NMI_inference_parameters( parsM= res$dstatM, parsrepM= res$dstatrepM , 
					fayfac=fayfac , RR=RR , Nimp=BIFIE.method$Nimp , 
					Nimp_NMI=BIFIE.method$Nimp_NMI , comp_cov = FALSE )									
		dfr$d <- res1$pars
		dfr$d_SE <- res1$pars_se
		dfr$d_t <- round( dfr$d / dfr$d_SE , 2 )
		dfr$d_df <- res1$d_df
		dfr$d_p <- stats::pt( - abs( dfr$d_t ) , df=dfr$d_df) * 2			
		dfr$d_fmi <- res1$pars_fmi
		dfr$d_VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
		dfr$d_VarRep <- res1$pars_varWithin	
						}		

	if ( ( ! res5$se ) &  ( RR==0 ) ){						
	# if (RR==0){				
		dfr$d_SE <- dfr$d_fmi <- dfr$d_VarMI <- dfr$d_VarRep <- NULL
				}

	if ( is_pseudogroup	){
			stat <- res5$stat[ 1:GG , ]
			stat$pseudogroup <- 1:GG
			ind1 <- grep( "groupvar" , colnames(stat) )
			ind2 <- grep( "groupval" , colnames(stat) )

			groupvar_pseudo <- apply( stat[ , ind1 ] , 1 , FUN = function(hh){ paste0( hh , collapse="#") } )
			groupval_pseudo <- apply( stat[ , ind2 ] , 1 , FUN = function(hh){ paste0( hh , collapse="#") } )
			dfr$group <- groupvar_pseudo[1]
			dfr$groupval1 <- groupval_pseudo[ dfr$groupval1 ]		
			dfr$groupval2 <- groupval_pseudo[ dfr$groupval2 ]
		
						}

			
	stat.dstat <- dfr	
	
	#*****
	# F statistics
	dfr <- NULL
	for (vv in 1:VV){
		Cdes <- matrix( 0 , nrow=GG-1 , ncol=GG*VV )
		indvec <- 1:(GG-1)
		for (zz in indvec ){
			Cdes[ zz , c(zz + (vv-1)*GG ,zz+1 + (vv-1)*GG ) ] <- c(1,-1)
					}
		rdes <- rep(0,GG-1)
		if ( wald_test ){
			wres5 <- BIFIE.waldtest( res5 , Cdes=Cdes , rdes=rdes )		
			dfr <- rbind( dfr , wres5$stat.D )
						}
				}
	if ( wald_test ){			
		dfr <- data.frame( "variable" = vars , "group" = group , dfr )	
		stat.F <- dfr
			} else {
		stat.F <- NA
				}
	
	
	parnames <- NULL
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat.F"= stat.F , "stat.eta" = stat.eta2 , "stat.dstat" = stat.dstat , 
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"NMI" = BIFIE.method$NMI , "Nimp_NMI" = BIFIE.method$Nimp_NMI , 
			"GG"=GG , "parnames" = parnames , "CALL"=cl , wald_test = wald_test )
	class(res1) <- "BIFIE.univar.test"
	return(res1)
		}
###################################################################################

####################################################################################
# summary
summary.BIFIE.univar.test <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	if ( object$wald_test ){
		cat("F Test (ANOVA) \n")	
		obji <- object$stat.F
		print.object.summary( obji , digits=digits )	
							}
	cat("\nEta Squared \n")	
	obji <- object$stat.eta
	print.object.summary( obji , digits=digits )
	cat("\nCohen's d Statistic \n")	
	obji <- object$stat.dstat
	print.object.summary( obji , digits=digits )	
			}