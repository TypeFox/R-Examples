

#######################################################################
# BIFIE.by function
BIFIE.by <- function( BIFIEobj , vars , userfct , userparnames=NULL ,
		group=NULL , group_values=NULL , se=TRUE , use_Rcpp = TRUE ){
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
						which( varnames == vv ) } , simplify=TRUE ) )
    # vars values
	VV <- length(vars)
					
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

		
	#****
	# pure R implementation
	if ( ! use_Rcpp ){	
		res <- BIFIE_by_helper_pureR(
			group_values , userfct , datalistM ,
			N , vars_index , wgt_ , wgtrep , Nimp , RR , fayfac ,
			group_index , userparnames
				)
					}
	
	#****
	# Rcpp implementation
	if ( use_Rcpp ){
		res <- .Call("bifie_by" , datalistM , wgt_ , wgtrep ,	vars_index - 1,    fayfac ,
				Nimp , group_index - 1 , group_values , userfct , PACKAGE="BIFIEsurvey")
					}
	
	NP <- res$NP
	GG <- length(group_values)
	ZZ <- NP
	if (is.null( userparnames ) ){
		userparnames <- paste0("parm",1:NP) 
				}
	
	dfr <- data.frame( "parm" = rep( userparnames , GG )
						)
	if (! nogroup){
	   dfr$groupvar <- group
	   dfr$groupval <- rep( group_values , each=ZZ )
	             }				 				 	


	dfr$Ncases <- rep( rowMeans( res$ncasesM ) , each=ZZ )
	dfr$Nweight <- rep( rowMeans( res$sumwgtM ) , each=ZZ )

	dfr <- create_summary_table( res_pars=res$parsL , 
				     parsM=res$parsM   , parsrepM=res$parsrepM , 
					 dfr=dfr , BIFIEobj=BIFIEobj )				
	dfr <- clean_summary_table( dfr=dfr , RR=RR , se=se , Nimp=Nimp )	
	
				
	# create vector of parameter names
	parnames <- paste0( dfr$parm   , "_" , dfr$groupvar , dfr$groupval )


	#@@@@***
	# multiple groupings
	dfr <- BIFIE_table_multiple_groupings( dfr , res00 )
	#@@@@***
						
	
	#*************************** OUTPUT ***************************************
	s2 <- Sys.time()
	timediff <- c( s1 , s2 ) # , paste(s2-s1 ) )
	res1 <- list( "stat" = dfr , 
			"output" = res , 	"timediff" = timediff ,
			"N" = N , "Nimp" = Nimp , "RR" = RR , "fayfac"=fayfac , "GG"=GG ,			
			"NMI" = BIFIEobj$NMI , "Nimp_NMI" = BIFIEobj$Nimp_NMI , 
			"parnames" = parnames , "CALL"= cl)
	class(res1) <- "BIFIE.by"
	return(res1)
		}
###################################################################################

####################################################################################
# summary for BIFIE.by function
summary.BIFIE.by <- function( object , digits=4 , ... ){
    BIFIE.summary(object)
	cat("Statistical Inference for User Defined Function \n")	
	obji <- object$stat
	print.object.summary( obji , digits=digits )
			}