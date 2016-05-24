#################################################
write.datlist <- function( datlist, name , include.varnames = TRUE ,
		type = "csv2" , separate = TRUE , Mplus = FALSE , round = NULL ,
		Rdata = TRUE , ... ){

		mplus <- Mplus
		
		# convert into datlist
		datlist <- mids2datlist(datlist)
		if ( inherits(datlist,"imputationList") ){
			datlist <- datlist$imputations
				}
		col.names <- include.varnames
		if (mplus){
			col.names <- FALSE
			type <- "table"
			}		
				
		# create subdirectory
		pf.subf <- file.path( getwd() , paste( name , sep=""))		
        dir.create(pf.subf)                 # define subdirectory

        # write legend of variables
		vars <- colnames(datlist[[1]])
		V <- length(vars)
		M <- length(datlist)
		dfr <- data.frame("index" = 1:V , "variable" = vars)
		
#        base::writeLines( vars , 
#		      file.path( pf.subf , paste( name , "__LEGEND.txt" , sep="") ))
        utils::write.table( dfr , 
				file.path( pf.subf , paste( name , "__LEGEND.txt" , sep="") ),
				row.names=FALSE , quote=FALSE
					)

		# declare types here!!	  
		type2 <- type	  
		if ( type=="csv2"){ type2 <- "csv"}
		if ( type=="table"){ type2 <- "dat"}		
		
        l1 <- paste( name , "__IMPDATA" , 1:M , "." , type2 , sep="")
        utils::write.table( l1 , file.path( pf.subf , 
		                  paste( name , "__IMP_LIST.txt" , sep="") ) , 
				       col.names=FALSE , row.names=FALSE , quote=FALSE)

		# save list of imputed datasets in Rdata format if requested
		if ( Rdata ){
			base::save( datlist , file=file.path( pf.subf ,
					paste( name , "__DATLIST.Rdata" , sep="") ) )
						}		
		# write datlist in separate files	
		if (separate){
		for (ii in 1:M){
			# ii <- 1
			file_ii <- paste( name , "__IMPDATA" , ii , sep="")
			dat_ii <-  datlist[[ii]]
			if ( ! is.null(round) ){
				dat_ii <- round( dat_ii , round)
					}
			args <- list( data = dat_ii , file = file_ii , path = pf.subf ,
							type = type , row.names= FALSE , 
							col.names=col.names, ... )
			if ( type %in% c("csv" , "csv2") ){
				args$col.names <- NULL
						}
			do.call( save.data , args )			
#			save.data( dat_ii , file = file_ii , path = pf.subf ,
#							type = type , row.names= FALSE , 
#							col.names=col.names, ... )	
			cat("Saved dataset" , ii ,"\n")
			utils::flush.console()
						}
				}

		#-----------
        # Mplus body                
        if ( mplus ){
		
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

				}
###################################################				