
######################################################################
# save BIFIEdata objects
save.BIFIEdata <- function( BIFIEdata , name.BIFIEdata , cdata=TRUE , varnames=NULL ){
    # bifieobj <- BIFIEdata	    
	make.cdata <- cdata
	if ( BIFIEdata$cdata ){ 
		make.cdata <- FALSE
			}			
    #********** convert BIFIE data object into compact BIFIEcdata object
	if ( make.cdata ){
		BIFIEdata <- BIFIE.BIFIEdata2BIFIEcdata( BIFIEdata , varnames = varnames )
		varnames <- NULL
				     }
	#********** compact saving of replicate weights	
	if (cdata){
		BIFIEdata$wgtreplist <- cdata.wgtrep( BIFIEdata$wgtrep )		
		BIFIEdata$wgtrep <- NULL
        BIFIEdata$cdata <- TRUE		
			  }	
    #******** variable selection in case of non-compact data
    if ( ! cdata ){	
			BIFIEdata <- BIFIE.data.select( BIFIEdata , varnames = varnames , 
						impdata.index = NULL )					
				}
	#******** variable selection in case of compact BIFIEdata
    if ( cdata ){		
		BIFIEdata <- BIFIE.cdata.select( BIFIEdata , 
					    varnames = varnames , impdata.index = NULL )							
				}		
	#***********************************************************
	#****** save objects
    base::save( BIFIEdata , file= paste0( name.BIFIEdata , ".Rdata") )
    sink( paste0( name.BIFIEdata , "__SUMMARY.Rout") )
		cat( getwd() , "\n" ,
			paste0( name.BIFIEdata , ".Rdata\n") , 
			"Saved at ", paste(Sys.time()) , 
			"\n\n")
		summary( BIFIEdata )
		cat("\n\nSaved variables:\n")
		VV <- length(BIFIEdata$varnames)
		for (vv in 1:VV){ 
			cat( vv , " " , BIFIEdata$varnames[vv] , "\n")
				}
    sink()
    cat( " - Saved" , paste0( name.BIFIEdata , ".Rdata") , "in directory \n    ")
    cat( getwd() , "\n")
        }
######################################################################
		
