##################################################
# write mice imputation object
write.mice.imputation <- function( mi.res , name , include.varnames = TRUE , long = TRUE , 
                           mids2spss = TRUE , spss.dec = "," , dattype = NULL ){

		
		ismids <- TRUE
		if ( inherits(mi.res ,"mids.1chain" ) ) {
                mi.res0 <- mi.res		
				mi.res <- mi.res$midsobj 
				mi.res$chainMean <- mi.res0$chainMean
				mi.res$chainVar <- mi.res0$chainVar					
				ismids <- FALSE
						}
		
        # pf.subf <- file.path( getwd() , paste("IMP_" , name , sep=""))
		pf.subf <- file.path( getwd() , paste( name , sep=""))
		
        dir.create(pf.subf)                 # define subdirectory
        # write legend of variables
        base::writeLines( colnames(mi.res$data) , 
		      file.path( pf.subf , paste( name , "__LEGEND.txt" , sep="") ))
        l1 <- paste( name , "__IMPDATA" , 1:mi.res$m , ".dat" , sep="")
        utils::write.table( l1 , file.path( pf.subf , 
		                  paste( name , "__IMP_LIST.txt" , sep="") ) , 
				       col.names=F , row.names=F , quote=F)
        # save summary in subdirectory
        sink( file.path( pf.subf , paste( name , "__IMP_SUMMARY.txt" , sep="") ) , split=TRUE )
        cat( paste(Sys.time()) , "\n\n" , pf.subf , "\n\n" )
            print( summary( mi.res ) ) 
        cat("\n\n") ; 
        ####
	    if ( ismids ){
        if ( ( mi.res$m > 1 ) & ( dim(mi.res$chainMean)[2]  > 1 )){    
				h1 <- Rhat.mice( mi.res ) 
				for (vv in seq(2,ncol(h1))){ h1[,vv] <- round( h1[,vv] , 3 ) }
				print(h1)
						} }
        print( utils::citation()) ; 
        print( utils::citation( "mice" ) ); 
        print(Sys.info()) ; 
        print( utils::sessionInfo())      
            sink()          ## end mice summary
        for (i in 1:mi.res$m ){ 
                        utils::write.table( mice::complete( mi.res , action = i ) , 
                                    file = file.path( pf.subf , paste( name , "__IMPDATA" , i , ".dat" , sep="")) , quote=F , 
                                        row.names=F , col.names=include.varnames , na ="." ) 
                if ( ! is.null(dattype) ){ 
                    if (dattype == "csv2" ){
                            write.csv2( mice::complete( mi.res , action = i ) , 
                                        file = file.path( pf.subf , paste( name , "__IMPDATA" , i , ".csv" , sep="")) , quote=F , 
                                            row.names=F ,  na ="." ) 
                                        }
                                    }
                        cat("\n",i); flush.console() 
                                }
        cat("\n")
        # write data file in long format
        utils::write.table( mice::complete( mi.res , action = "long" ) , 
                                    file = file.path( pf.subf , paste( name , "__LONG.dat" , sep="")) , quote=F , 
                                        row.names=F , col.names= TRUE , na ="." )         
        # variable names
        base::writeLines( colnames( mice::complete( mi.res , 1 ) ) , file.path( pf.subf , paste( name , "__VARNAMES.txt" , sep="")) )
        # write SPSS file
        if (mids2spss){ 
            mice::mids2spss(mi.res, filedat = paste( name , "__SPSS.txt" , sep="") , 
                               filesps = paste( name , "__SPSS.sps", sep="") , 
                               path = pf.subf , dec = spss.dec )
                                }
        # Mplus-Body                
        if ( ! include.varnames ){
		
			vars <- colnames(mi.res$data) 		
			vars2 <- VariableNames2String( vars , breaks=60)
		
            l1 <- c("TITLE: xxxx ;" , "" , "DATA: " , "" ,
                    paste( "FILE IS " , name , "__IMP_LIST.txt;" , sep=""),  "TYPE = IMPUTATION;"  , "" , 
                        "VARIABLE:" , "" , "NAMES ARE" , 
                       vars2 , ";" , "" , "! edit usevariables are;" , "!usevar are" , 
                        "   " , "" , "MISSING = ." , "" , "!........................." ,
                        "! Mplus statements"
                            )
            base::writeLines( l1 , file.path( pf.subf , paste( name , "__MPLUS-INPUT-BODY.inp" , sep="") ))
            }
        # write predictorMatrix and imputationMethod
        utils::write.csv2( mi.res$method  , file.path( pf.subf , paste( name , "__IMPMETHOD.csv" , sep="")) , quote=F)
        utils::write.csv2( mi.res$predictorMatrix  , file.path( pf.subf , paste( name , "__PREDICTORMATRIX.csv" , sep="")) , quote=F)               
		# write mice object as a Rdata object
		base::save( mi.res , file=file.path( pf.subf , paste( name , ".Rdata" , sep="") ) )
		# save list of imputed datasets
		datlist <- mids2datlist( mi.res )
		base::save( datlist , file=file.path( pf.subf , paste( name , "__DATALIST.Rdata" , sep="") ) )		
        }
#########################################################################