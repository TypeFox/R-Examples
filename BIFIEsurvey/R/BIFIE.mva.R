

#######################################################################
# Missing value analysis
BIFIE.mva <- function( BIFIEobj , missvars , covariates=NULL , se=TRUE ){
	#****
	s1 <- Sys.time()
	bifieobj <- BIFIEobj
	cl <- match.call()		
#	if (bifieobj$cdata){
#		varnames <- unique( c( vars , group , "one") )
#		bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames=varnames )	
#						}	
	# 
	if ( ! bifieobj$cdata ){
		varnames <- unique( c(missvars , covariates ) )
		bifieobj <- BIFIE.BIFIEdata2BIFIEcdata( bifieobj , varnames )
						}

	if ( is.null(covariates) ){
		covariates <- "one"
					}
						
	# if covariates = NULL, then create a garbage variable
	if ( is.null(covariates) ){ 
#	if ( FALSE ){ 
		N <- bifieobj$N
		transform.formula <-  paste0( "~ 0 + I ( stats::runif( " , N , " , 0 , 1E-10) ) " )
#		transform.formula <-  paste0( "~ 0 + I ( 1: ", N , "  ) " )
		bifieobj <- BIFIE.data.transform( bifieobj , transform.formula ,  "_null" )  
		covariates <- bifieobj$varnames.added
		se <- FALSE
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
		
	# start with a compact BIFIEdata object
	varnames <- unique( c(missvars , covariates ) )
#	varnames <- setdiff( varnames , "one" )

	# assume that bifieobj is already in cdata format
	
	# define selected response indicators

	missvars_index <- match( missvars , colnames(bifieobj$dat1) )
	datalistM_ind_sel <- bifieobj$datalistM_ind[ , missvars_index , drop=FALSE ]											

	respvars <- paste0("resp_" , missvars )
	colnames(datalistM_ind_sel) <- respvars
	VVadd <- length(respvars )

	varnames1 <- c( varnames , respvars )

	bifieobj$datalistM_ind <- cbind( bifieobj$datalistM_ind , datalistM_ind_sel )
	bifieobj$dat1 <- cbind( bifieobj$dat1 , datalistM_ind_sel )
	bifieobj$Nvars <- bifieobj$Nvars + VVadd
	bifieobj$varnames <- c( bifieobj$varnames , respvars )


	# select dataset
	bifieobj <- BIFIE.BIFIEcdata2BIFIEdata( bifieobj , varnames = varnames1 )
	
	RR <- bifieobj$RR
    do_test <- TRUE
    if (RR < 2 ){
		wgtrep <- bifieobj$wgtrep 
		do_test <- FALSE
		wgtrep <- cbind( wgtrep , wgtrep + stats::runif( nrow(wgtrep) , 0 , 1E-4) )
		bifieobj$wgtrep <- wgtrep
		bifieobj$RR <- 2
				}
	
	if ( is.null( covariates) ){
			do_test <- FALSE 
			se <- FALSE
			covariates <- "one"
				}
	
	if (RR==1){ RR <- 0 }
	if ( ! se ){ 
		N11 <- length(wgt)
		wgtrep <- matrix( NA , nrow=N11 , ncol=2 )
		wgtrep[,1] <-wgt
		eps <- 1E-8
		wgtrep[,2] <- wgt + stats::runif(N11, -eps,eps)
		RR <- 0
		bifieobj$wgtrep <- wgtrep
		bifieobj$RR <- ncol(wgtrep)
				}	
		
	#*****
	# collect results
	VV <- length(respvars)
	CVV <- length(covariates)
	res_list <- list( 1:VV)
	dfr <- NULL
	
	RR <- bifieobj$RR

	
	for (vv in 1:VV){
		# vv <- 1
		rvv <- respvars[vv]
		res.vv <- as.list(1:2)
		names(res.vv) <- c("stat" , "dstat")	
		res <- BIFIE.univar( bifieobj , vars = covariates , group= rvv, 
					group_values = 0:1 )
		res.vv$stat <- res$stat
		
		res1 <- BIFIE.univar.test( res , wald_test = FALSE )					
		res.vv$dstat <- res1$stat.dstat
						
		res_list[[vv]] <- res.vv
		# collect results
		dfr.vv <- data.frame( "respvar"= rep(rvv,CVV) )		
		dfr.vv$missprop <- res$stat$Nweight[1] / ( res$stat$Nweight[1] + res$stat$Nweight[2] )
		dfr.vv$covariate <- covariates
		if ( do_test ){
			dfr.vv$d <- res1$stat.dstat$d
			dfr.vv$d_SE <- res1$stat.dstat$d_SE
			dfr.vv$t <- res1$stat.dstat$t
			dfr.vv$p <- res1$stat.dstat$p
						}						
		dfr.vv$M_resp <- res$stat$M[ seq(2,2*CVV , 2 ) ]
		dfr.vv$M_miss <- res$stat$M[ seq(1,2*CVV , 2 ) ]
		dfr.vv$SD_resp <- res$stat$SD[ seq(2,2*CVV , 2 ) ]
		dfr.vv$SD_miss <- res$stat$SD[ seq(1,2*CVV , 2 ) ]
		dfr <- rbind( dfr , dfr.vv )
				}
	if ( covariates[1] == "_null" ){ se <- FALSE }
	
	if ( ! do_test ){
		RR <- 0
				  }
	
	if ( ( ! se ) |  ( RR==0 ) ){
		dfr$t <- dfr$p <- dfr$d_SE <-  NULL
				}						
	
	if ( covariates[1] %in% c("_null","one") ){			
			dfr$covariate <- dfr$d <- dfr$M_resp  <- NULL
			dfr$M_miss <- dfr$SD_resp <- dfr$SD_miss <- NULL						
			}
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat.mva" = dfr , "res_list" = res_list , 
			"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac ,
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
			"CALL"= cl
			# "itempair_index" = itempair_index , "GG"=GG ,
			# "parnames" = parnames
				)
	class(res1) <- "BIFIE.mva"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.mva function
summary.BIFIE.mva <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Missing Value Analysis \n")	
	obji <- object$stat.mva
	print.object.summary( obji , digits=digits )			
			}