mids2datlist <- function( midsobj , X = NULL){
    datlist <- midsobj   # init
	if ( class(midsobj) == "mids.1chain" ){
		midsobj <- midsobj$midsobj		
			}
	#*****************************		
	# object of class mids
	if ( class(midsobj) %in% c("mids" ) ){		
		m <- midsobj$m
		datlist <- as.list( 1:m )
		for (ii in 1:m){  
				h1 <- mice::complete( midsobj , ii ) 
				if ( ! is.null(X) ){
					h1 <- data.frame(  X , h1	)										
								   }	
				datlist[[ii]] <- h1
						}
		class(datlist) <- "datlist"
		datlist <- datlist_create(datlist)
					}	
	#******************************
    # object of class mids.nmi		
	if ( class(midsobj) %in% c("mids.nmi" ) ){		
		Nimp <- midsobj$Nimp
		datlist <- as.list(1:Nimp["between"])
		dat1 <- as.list(1:Nimp["within"])
		# imp <- midsobj$imp
		for (bb in 1:Nimp["between"]){
			# bb <- 1
			datlist[[bb]] <- dat1
			for (ww in 1:Nimp["within"] ){
				# ww <- 1
				h1 <- complete.mids.nmi( midsobj , action= c(bb,ww) )	
				if ( ! is.null(X) ){
					h1 <- data.frame(  X , h1	)										
								   }				
				datlist[[bb]][[ww]] <- h1
										}
								}
		class(datlist) <- "nested.datlist"
		datlist <- nested.datlist_create(datlist)
			}		
    return(datlist)
        }
