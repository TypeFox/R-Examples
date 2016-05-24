
#######################################################################
# univariate statistics
BIFIE.univar <- function( BIFIEobj , vars , group=NULL , group_values=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	cl <- match.call()		
	bifieobj <- BIFIEobj
	if (bifieobj$cdata){
		varnames <- unique( c( vars , group , "one") )
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
	
	vars_index <- unlist( sapply( vars , FUN = function(vv){ 
						which( varnames == vv ) } ) )		
	wgt_ <- matrix( wgt , ncol=1 )
	if ( is.null( group) ){ 
			nogroup <- TRUE } else { 
			nogroup <- FALSE 
#			if (group=="one"){ nogroup <- TRUE }
					}	
					
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
				
	#**************************************************************************#
	#****************** no grouping variable **********************************#
	if ( nogroup ){
#		res <- .Call( "univar_1group" ,  datalistM , wgt_ , wgtrep , vars_index - 1 , fayfac ,  Nimp ,
#				PACKAGE="BIFIEsurvey" ) 
#		GG <- 1
#		VV <- length(vars)

		res <- .Call( "univar_multiple_V2group" , datalistM , wgt_ , wgtrep , vars_index - 1 , fayfac , Nimp ,
				group_index - 1, group_values , PACKAGE="BIFIEsurvey" )								
		GG <- length(group_values)
		VV <- length(vars)					
		dfr <- data.frame( "var" = rep(vars,each=GG) , 
				"Nweight" = rowMeans(res$sumweightM) ,
				"Ncases" = rowMeans( res$ncasesM) , 
				"M" = res$mean1 , "M_SE" = res$mean1_se )		
		dfr$M_df <- round( (Nimp-1)*( 1 +  (Nimp*res$mean1_varWithin )/ ( Nimp+1) / res$mean1_varBetween )^2  , 2 )
		vv <- "M_df"
		dfr[,vv] <- ifelse( dfr[,vv] > 1000 , Inf , dfr[,vv] )
		
		dfr$M_t <- dfr$M / dfr$M_SE
		dfr$M_p <- 2* stats::pt( - abs( dfr$M_t) , df = dfr$M_df )			
		dfr0 <- data.frame(		"M_fmi" = res$mean1_fmi   , 
		        "M_VarMI"= res$mean1_varBetween , "M_VarRep"= res$mean1_varWithin   ,
				"SD" = res$sd1 , "SD_SE" = res$sd1_se )	
		dfr <- cbind( dfr ,  dfr0 )				
		dfr$SD_df <- round( (Nimp-1)*( 1 +  (Nimp*res$sd1_varWithin )/ ( Nimp+1) / res$sd1_varBetween )^2  , 2 )
		vv <- "SD_df"
		dfr[,vv] <- ifelse( dfr[,vv] > 1000 , Inf , dfr[,vv] )

		dfr$SD_t <- dfr$M / dfr$SD_SE
		dfr$SD_p <- 2*stats::pt( - abs( dfr$SD_t) , df = dfr$SD_df )							
		dfr0 <- data.frame(		"SD_fmi" = res$sd1_fmi   , "SD_VarMI"= res$sd1_varBetween , "SD_VarRep"= res$sd1_varWithin            
					)

				
		if (BIFIEobj$NMI ){
			# M
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$mean1M , parsrepM=res$mean1repM , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )						
			dfr$M <- res1$pars
			dfr$M_SE <- res1$pars_se
			dfr$M_df <- res1$df
			dfr$M_t <- res1$pars / res1$pars_se 
			dfr$M_p <- 2*stats::pt( - abs( dfr$M_t) , df = res1$df )			
			dfr$M_fmi <- res1$pars_fmi
			dfr$M_VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr$M_VarRep <- res1$pars_varWithin		
			# SD			
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$sd1M , parsrepM=res$sd1repM , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )			
			dfr$SD <- res1$pars
			dfr$SD_SE <- res1$pars_se
			dfr$SD_df <- res1$df
			dfr$SD_t <- res1$pars / res1$pars_se 
			dfr$SD_p <- 2*stats::pt( - abs( dfr$SD_t) , df = res1$df )						
			dfr$SD_fmi <- res1$pars_fmi
			dfr$SD_VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr$SD_VarRep <- res1$pars_varWithin									
							}
						
				
				}
	#**************************************************************************#
	#****************** with grouping variable ********************************#		
	if ( ! nogroup ){

		res <- .Call( "univar_multiple_V2group" , datalistM , wgt_ , wgtrep , vars_index - 1 , fayfac , Nimp ,
				group_index - 1, group_values , PACKAGE="BIFIEsurvey" )	
		GG <- length(group_values)
		VV <- length(vars)
					 
		dfr <- data.frame( "var" = rep(vars,each=GG) , 
				"groupvar" = group , 
				"groupval" = rep(group_values  , VV ) , 
				"Nweight" = rep( rowMeans(res$sumweightM) , VV ) ,
				"Ncases" = res$ncases ,
				"M" = res$mean1 , "M_SE" = res$mean1_se )
		dfr$M_df <- round( (Nimp-1)*( 1 +  (Nimp*res$mean1_varWithin )/ ( Nimp+1) / res$mean1_varBetween )^2  , 2 )
		vv <- "M_df"
		dfr[,vv] <- ifelse( dfr[,vv] > 1000 , Inf , dfr[,vv] )

		dfr$M_t <- dfr$M / dfr$M_SE
		dfr$M_p <- 2*stats::pt( - abs( dfr$M_t) , df = dfr$M_df )			
		dfr <- data.frame( dfr ,	"M_fmi" = res$mean1_fmi   , "M_VarMI"= res$mean1_varBetween , "M_VarRep"= res$mean1_varWithin   ,
				"SD" = res$sd1 , "SD_SE" = res$sd1_se )
		dfr$SD_df <- round( (Nimp-1)*( 1 +  (Nimp*res$sd1_varWithin )/ ( Nimp+1) / res$sd1_varBetween )^2  , 2 )
		vv <- "SD_df"
		dfr[,vv] <- ifelse( dfr[,vv] > 1000 , Inf , dfr[,vv] )
		dfr$SD_t <- dfr$M / dfr$SD_SE
		dfr$SD_p <- 2*stats::pt( - abs( dfr$SD_t) , df = dfr$SD_df )							
		dfr <- data.frame( dfr  ,	"SD_fmi" = res$sd1_fmi   , "SD_VarMI"= res$sd1_varBetween , "SD_VarRep"= res$sd1_varWithin            
					 )
						
		if (BIFIEobj$NMI ){
			# M
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$mean1M , parsrepM=res$mean1repM , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )			
			dfr$M <- res1$pars
			dfr$M_SE <- res1$pars_se
			dfr$M_df <- res1$df
			dfr$M_t <- res1$pars / res1$pars_se 
			dfr$M_p <- 2*stats::pt( - abs( dfr$M_t) , df = res1$df )						
			dfr$M_fmi <- res1$pars_fmi
			dfr$M_VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr$M_VarRep <- res1$pars_varWithin		
			# SD			
			res1 <- BIFIE_NMI_inference_parameters( parsM=res$sd1M , parsrepM=res$sd1repM , 
						fayfac=fayfac , RR=RR , Nimp=Nimp , 
						Nimp_NMI=BIFIEobj$Nimp_NMI , comp_cov = FALSE )			
			dfr$SD <- res1$pars
			dfr$SD_SE <- res1$pars_se
			dfr$SD_df <- res1$df
			dfr$SD_t <- res1$pars / res1$pars_se 
			dfr$SD_p <- 2*pt( - abs( dfr$SD_t) , df = res1$df )				
			dfr$SD_fmi <- res1$pars_fmi
			dfr$SD_VarMI <- res1$pars_varBetween1 + res1$pars_varBetween2
			dfr$SD_VarRep <- res1$pars_varWithin									
							}
					 					 
	         	}
								
				
	if ( ( ! se ) &  ( RR==0 ) ){				
		dfr$M_SE <- dfr$M_fmi <- dfr$M_VarMI <- dfr$M_VarRep <- dfr$M_t <- dfr$M_df <- dfr$M_p <- NULL
		dfr$SD_SE <- dfr$SD_fmi <- dfr$SD_VarMI <- dfr$SD_VarRep <- 
					dfr$SD_t <- dfr$SD_df <- dfr$SD_p <-NULL		
				}	
	if ( Nimp == 1 ){				
		dfr$M_fmi <- dfr$M_VarMI <- NULL
		dfr$SD_fmi <- dfr$SD_VarMI <- NULL		
				}		
				
	
	#****
	# statistics for mean and SD
	v1 <- c("var" , "groupvar" , "groupval" , "Nweight" , "Ncases" )
	v1 <- intersect( v1 , colnames(dfr) )
	cdfr <- colnames(dfr)
	stat_M <- dfr[ , c( v1  , cdfr[ substring( cdfr , 1 ,1) == "M" ] ) ]
	stat_SD <- dfr[ , c( v1  , cdfr[ substring( cdfr , 1 ,2) == "SD" ] ) ]
	
	# create vector of parameter names
	nogroupL <- rep( nogroup , nrow(dfr) )
	parnames <- paste0( dfr$var   , 
			ifelse( ! nogroupL , paste0( "_" , dfr$groupvar , "_" ) , "" ) ,
			ifelse( ! nogroupL , dfr$groupval , "" ) )		

	#@@@@***
	# multiple groupings
	dfr <- BIFIE_table_multiple_groupings( dfr , res00 )
	#@@@@***
	
	stat_M <- BIFIE_table_multiple_groupings( stat_M , res00 )
	stat_SD <- BIFIE_table_multiple_groupings( stat_SD , res00 )	
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr , "stat_M" = stat_M , "stat_SD" = stat_SD , 
			"output" = res , "timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac , "parnames" = parnames ,
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , "se"=se ,
			"GG" = GG , "VV"=VV , "vars" = vars , "group" = group , "CALL"=cl)
	class(res1) <- "BIFIE.univar"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.univar function

summary.BIFIE.univar <- function( object , digits=3 , ... ){
    BIFIE.summary(object)     
	cat("Univariate Statistics | Means\n")	
	obji <- object$stat_M
	print.object.summary( obji , digits=digits )
	cat("\nUnivariate Statistics | Standard Deviations\n")	
	obji <- object$stat_SD
	print.object.summary( obji , digits=digits )
	
			}