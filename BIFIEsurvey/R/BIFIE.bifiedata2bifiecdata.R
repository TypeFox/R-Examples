
#######################################################
# conversion of BIFIEdata to BIFIEcdata
BIFIE.BIFIEdata2BIFIEcdata <- function( bifieobj , varnames=NULL , impdata.index = NULL ){					
   if ( bifieobj$cdata ){
		stop( "You may want to use 'BIFIE.BIFIEcdata2BIFIEdata'\n")
						}
	#******** select some imputed datasets or some variables
	bifieobj <- BIFIE.data.select( bifieobj=bifieobj , varnames=varnames , 
					       impdata.index =impdata.index )
						   
	#**** data conversion
    res1 <- .Call("bifie_bifiedata2bifiecdata" ,
					    bifieobj$datalistM , bifieobj$Nimp ,
					    PACKAGE="BIFIEsurvey" )
    bifieobj$cdata <- TRUE
    bifieobj$datalistM <- NULL
    bifieobj$datalistM_ind <- res1$datalistM_ind	
	colnames(bifieobj$datalistM_ind) <- bifieobj$varnames
    bifieobj$datalistM_imputed <- res1$datalistM_imputed
	bifieobj$datalistM_impindex <- res1$datalistM_impindex	
#	colnames(bifieobj$datalistM_imputed) <- c("_imp" , "subj" , "variable" , "value")
    bifieobj$time <- Sys.time()
    return(bifieobj)
            }
#######################################################			

#######################################################
# conversion of BIFIEcdata to BIFIEdata object
BIFIE.BIFIEcdata2BIFIEdata <- function( bifieobj , varnames=NULL , impdata.index = NULL ){
   if ( ! bifieobj$cdata ){
		stop( "You may want to use 'BIFIE.BIFIEdata2BIFIEcdata'\n")
						}											
	#******** select some imputed datasets or some variables
	bifieobj <- BIFIE.cdata.select( bifieobj=bifieobj , varnames=varnames , 
					impdata.index =impdata.index )				

	#***** conversion to BIFIEdata object
	bifieobj$datalistM <- .Call("bifie_bifiecdata2bifiedata" ,
								   as.matrix(bifieobj$datalistM_ind) , 
								   as.matrix(bifieobj$datalistM_imputed) , 
								   bifieobj$Nimp , 
								   as.matrix(bifieobj$dat1) ,
								   as.matrix(bifieobj$datalistM_impindex) , 
								   PACKAGE="BIFIEsurvey" )$datalistM
	bifieobj$cdata <- FALSE
	bifieobj$datalistM_imputed <- NULL
    bifieobj$datalistM_impindex <- NULL	
	bifieobj$datalistM_ind <- NULL	
	return(bifieobj)
}
##############################################################