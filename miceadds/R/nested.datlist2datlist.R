
##########################################################
# converts a nested list of multiply imputed
# datasets into a list of multiply imputed datasets
nested.datlist2datlist <- function(datlist){
		CALL <- match.call()
		if (inherits(datlist,"NestedImputationList")){
			is_NestedImputationList <- TRUE
		    datlist <- datlist$imputations 
				} else {
			is_NestedImputationList <- FALSE					
				}
		Nimp <- c( length(datlist) , length(datlist[[1]] ) )
		names(Nimp) <- c("Between" , "Within")
		PP <- Nimp[1] * Nimp[2]
		datlist0 <- as.list( 1:PP )
		vv <- 1
		for (bb in 1:Nimp[1] ){
			for (ww in 1:Nimp[2] ){
					datlist0[[vv]] <- datlist[[bb]][[ww]]
					vv <- vv + 1
								}
							}
		# attr(datlist0,"Nimp") <- Nimp
		# class(datlist0) <- "datlist"
		if (! is_NestedImputationList ){
			datlist0 <- datlist_create(datlist0)
			attr(datlist0,"call") <- CALL
							} else {
			datlist0 <- mitools::imputationList(datlist0)
			datlist0$call <- CALL
								}							
		return(datlist0)						
				}
##########################################################				

####################################################
# datlist -> nested.datlist
datlist2nested.datlist <- function(datlist, Nimp){
		CALL <- match.call()
		if (inherits(datlist,"imputationList")){
			is_imputationList <- TRUE
		    datlist <- datlist$imputations 
				} else {
			is_imputationList <- FALSE					
				}

		PP <- Nimp[1] * Nimp[2]
		datlist1 <- as.list( 1:Nimp[1] )
		datlist2 <- as.list( 1:Nimp[2] )		
		vv <- 1
		for (bb in 1:Nimp[1] ){		
			for (ww in 1:Nimp[2] ){
					datlist2[[ww]] <- datlist[[vv]]
					vv <- vv + 1
								}
			datlist1[[bb]] <- datlist2				
							}
		if ( ! is_imputationList){					
			datlist1 <- nested.datlist_create(datlist1)					
			attr(datlist1,"call") <- CALL
								} else {								
			datlist1 <- NestedImputationList(datlist1)					
			datlist1$call <- CALL															
								}							
		return(datlist1)						
				}
#####################################################				