
##################################################################
##################################################################
##################################################################
# loading files for conversion into BIFIE data objects
load.BIFIEdata.files <- function( files.imp , wgt , file.wgtrep , 
        file.ind = NULL , type="Rdata" , 
		varnames=NULL , cdata =TRUE , dir=getwd() , ... ){
	Nimp <- length(files.imp)		
	# handle cases in which no weights are defined?
	# wgt <- "TOTWGT"			
	#########################################################
	# no indicator dataset
	#########################################################
	if ( is.null(file.ind) ){
		#**********
		# read imputed datasets
		datalist <- list(1:Nimp)
		for (ii in 1:Nimp){
			# ii <- 1
			cat(paste0( "- Read " , files.imp[[ii]] , "\n") ); utils::flush.console() ;
			dat1 <- miceadds::load.data( file = files.imp[[ii]] , path = dir , type=type , ... )
			if (ii==1){  wgt <- dat1[ , wgt ] }
			dat1 <- as.data.frame(dat1)
			if ( ! is.null(varnames) ){	
					dat1 <- dat1[ , varnames ]
						}		
			dat1$one <- NULL
			datalist[[ii]] <- dat1		
					}    		
		#**************
		# read replicate weights	
		cat(paste0( "- Read " , file.wgtrep , "\n") ); utils::flush.console() ;
		wgtrep <- miceadds::load.data( file = file.wgtrep , type=type , path = dir , ...)
		#****************************
		# create BIFIEdata object
		bifieobj <- BIFIE.data( data.list= datalist , wgt = wgt , wgtrep=wgtrep , 
                 						cdata=cdata )		
					}	
	#############################################################
	# with indicator dataset	
	#############################################################					
	if ( ! is.null( file.ind ) ){
		#*********************
		# Read indicator dataset
		cat(paste0( "- Read " , file.ind , "\n") ); utils::flush.console() ;
		dat_ind <- miceadds::load.data( file = file.ind , type=type , path = dir , ...)
		if ( is.null(varnames) ){
		     varnames <- setdiff( colnames(dat_ind ) , "one" )
								}
		if ( ! is.null(varnames) ){
			dat_ind <- dat_ind[ , varnames ]			
							}
		dat_ind <- as.matrix( dat_ind )
		# add column 1 for "one"
		dat_ind <- cbind( dat_ind , 1 )
		colnames(dat_ind) <- c( varnames , "one")
		
		#************************
		# Read first imputed dataset
		ii <- 1
		cat(paste0( "- Read " , files.imp[[ii]] , "\n") ); utils::flush.console() ;
		dat1 <- miceadds::load.data( file = files.imp[[ii]] , path = dir , type=type , ... )
		if ( ! is.null(varnames) ){
			dat1 <- dat1[ , varnames ]			
							}		
        datalist <- list( dat1 )
		#**************
		# read replicate weights	
		cat(paste0( "- Read " , file.wgtrep , "\n") ); utils::flush.console() ;
		wgtrep <- miceadds::load.data( file = file.wgtrep , type=type , path = dir , ...)
		#***************
		# create initial BIFIEdata object		
		bifieobj <- BIFIE.data( data.list= datalist , wgt = wgt , wgtrep=wgtrep , cdata=cdata )
        bifieobj$datalistM_ind <-  dat_ind	
        bifieobj$Nimp <- Nimp        		
		Nmiss <- sum( 1 - dat_ind )
		datalistM_imputed <- matrix( NA , nrow=Nmiss , Nimp)		
		# datalistM_impindexL <- NULL
		res1 <- .Call("bifie_bifiedata_stepwise" ,  as.matrix(dat1) , dat_ind , Nmiss ,
							PACKAGE="BIFIEsurvey" )$datalistM_imputed
		# bifieobj$datalistM_imputed <- res1		
		datalistM_imputed[,1] <- res1[,4]
		datalistM_impindex <- res1[,2:3]
		
		#************************
		# Read other imputed datasets
		if (Nimp>1){
		for (ii in 2:Nimp){
			cat(paste0( "- Read " , files.imp[[ii]] , "\n") ); utils::flush.console() ;
			dat1 <- miceadds::load.data( file = files.imp[[ii]] , path = dir , type=type , ... )
			if ( ! is.null(varnames) ){
				dat1 <- dat1[ , varnames ]			
								}		
			res1 <- .Call("bifie_bifiedata_stepwise" ,  as.matrix(dat1) , dat_ind , Nmiss ,
							PACKAGE="BIFIEsurvey" )$datalistM_imputed
#			res1[,1] <- ii-1
			datalistM_imputed[,ii] <- res1[,4]
		    datalistM_impindex <- rbind( datalistM_impindex , res1[,2:3] )
#			bifieobj$datalistM_imputed <- as.matrix( rbind( bifieobj$datalistM_imputed , res1 ) )
							}						
						}
        bifieobj$dat1 <- cbind( as.data.frame(dat1) , "one"=1 )
		bifieobj$datalistM_imputed <- datalistM_imputed	
		datalistM_imputed <- NULL
		bifieobj$datalistM_impindex <- datalistM_impindex		
#	    colnames(bifieobj$datalistM_imputed) <- c("_imp" , "subj" , "variable" , "value")		
					}  # end indicator data
    #####################################			
	# return BIFIE object
	return(bifieobj)	
		}
##################################################################
##################################################################
##################################################################
		