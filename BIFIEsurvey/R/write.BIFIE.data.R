
#############################################################
# write BIFIEdata object
write.BIFIEdata <- function( BIFIEdata , name.BIFIEdata , 
		dir = getwd() , varnames = NULL ,
		impdata.index = NULL , type = "Rdata" , ... ){	
	#****************************************
	dir1 <- getwd()	
	setwd(dir)
#	dir -> path.BIFIEdata
#    BIFIEdata <- BIFIEdata	
	cdata <- BIFIEdata$cdata			
	# make BIFIE object smaller
	BIFIEdata <- BIFIEdata.select( BIFIEdata , varnames  , impdata.index)
	cat("** Working directory:" , dir , "\n")
	#*************************
	# define file suffixes
	filesuf <- paste0("." , type )
	if ( type == "csv2" ){ filesuf <- ".csv" }
	if ( type == "table" ){ filesuf <- ".dat" }
	if ( type == "sav" ){ filesuf <- "" }
	#***************************************
	# save dataset with replicate weights
	cat(" - Saved replicate weights\n") ; utils::flush.console()	
	filename.temp <- paste0( name.BIFIEdata , "__WGTREP" , filesuf )
	w1 <- as.data.frame( BIFIEdata$wgtrep )
	miceadds::save.data( w1 , filename=filename.temp , type=type , path=dir , ... ) 
	w1 <- NULL							
	#*******************************************
	# save imputed datasets
	Nimp <- BIFIEdata$Nimp
	for (ii in 1:Nimp){
		# ii <- 1  
		cat(" - Saved imputed dataset" , ii , "\n") ; utils::flush.console()
		if (! cdata ){
			bii <- BIFIEdata.select( BIFIEdata , varnames = varnames , impdata.index = ii )
						}
		if ( cdata ){
			bii <- BIFIE.BIFIEcdata2BIFIEdata( BIFIEdata , varnames = varnames ,
			                impdata.index = ii )
						}						
		dat1 <- bii$datalistM
		colnames(dat1) <- bii$varnames
		filename.temp <- paste0( name.BIFIEdata , "__IMP" , ii , filesuf )
		dat1 <- as.data.frame(dat1)
	    miceadds::save.data( dat1 , filename=filename.temp , type=type , path=dir , ... ) 
			}
	#**********************************************
	# save BIFIEdata object
	save.BIFIEdata( BIFIEdata , paste0( name.BIFIEdata , "__BIFIEdataObject" ) , cdata = TRUE )	
	#**********************************************
	# finish
	setwd(dir1)	
				}