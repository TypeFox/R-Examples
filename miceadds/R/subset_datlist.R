
#######################################################
subset_datlist <- function( datlist , subset = TRUE ,
					select = NULL , expr_subset = NULL , index = NULL , 
					toclass = "datlist"){
		CALL <- match.call()
		#-----------------------------
		#*** check here for classes
		if ( class(datlist) %in% "imputationList" ){
			datlist <- datlist$imputations		
						}
		if ( class(datlist) %in% c("mids","mids.1chain") ){
			datlist <- mids2datlist( datlist )
									}		
		
		#*****************************
		# check for expr
		expr <- expr_subset 
		is_expr <- FALSE
		pf <- parent.frame()
		apply_select0 <- FALSE
	    if (!is.null(match.call()$expr)){
		    expr1 <- substitute(expr)
			is_expr <- TRUE
			apply_select0 <- TRUE
					}
																										
		#-----------------------------
		#*** start routine
		M <- length(datlist)				
		if ( is.null(index) ){
			index <- 1:M
						}		
		IM <- length(index)				
		if( is.null(select) & ( mean( subset ) == 1 ) ){
			apply_select <- FALSE
						} else {
			apply_select <- TRUE			
						}
		if (apply_select0){ apply_select <- TRUE }
		if ( is.null(select) ){
				select <- colnames(datlist[[1]])			
						}		
		datlist2 <- as.list(1:IM)
		for (ii in 1:IM){
			d1 <- datlist[[ index[ii] ]]	
			if (is_expr){
				 # h1 <- with( data=d1 , { expr1 } )
				 subset <- eval(expr1, d1, enclos=pf)
						}
			if (apply_select){

					d1 <- subset( d1 , subset=subset , select=select , drop=FALSE)
							 }
			datlist2[[ii]] <- d1
						}
		#************
		# create object classes
		#---- class datlist
		if (toclass == "datlist" ){
			datlist2 <- datlist_create(datlist2)
			attr(datlist2,"call") <- CALL
								   }	
		#---- class imputationList
		if (toclass == "imputationList" ){
			datlist2 <- mitools::imputationList(datlist2)
			datlist2$call <- CALL
								   }	
		#---- class mids
		if (toclass == "mids" ){
			datlist2 <- datlist2mids(datlist2, progress=TRUE)
			datlist2$call <- CALL
								}		
		return(datlist2)
		}
############################################################
# object of class datlist
subset.datlist <- function( x , subset ,
					select = NULL ,  expr_subset = NULL , 
					index = NULL , ... ){
		CALL <- match.call()
		if (missing(subset)){  subset <- TRUE	}
		datlist2 <- subset_datlist( datlist = x , subset = subset ,
					   select = select ,  expr_subset=expr_subset,
					   index = index , toclass = "datlist")					   
		attr(datlist2,"call") <- CALL
		return(datlist2)
				}
#---------------------------------------------------------------				
# object of class mids
subset.mids <- function( x , subset ,
					select = NULL ,  expr_subset = NULL ,  
					index = NULL , ... ){
		CALL <- match.call()
		if (missing(subset)){  subset <- TRUE }		
		datlist2 <- subset_datlist( datlist = x , subset = subset ,
					   select = select ,  expr_subset=expr_subset,
					   index = index , toclass = "mids")					   
		datlist2$call <- CALL
		return(datlist2)
				}
#---------------------------------------------------------------				
# object of class mids.1chain
subset.mids.1chain <- function( x , subset ,
					select = NULL ,  expr_subset = NULL , 
					index = NULL , ... ){
		CALL <- match.call()
		if (missing(subset)){  subset <- TRUE }		
		datlist2 <- subset_datlist( datlist = x , subset = subset ,
					   select = select ,  expr_subset=expr_subset,
					   index = index , toclass = "mids")					   
		datlist2$call <- CALL
		return(datlist2)
				}				
#---------------------------------------------------------------						
# object of class imputationList
subset.imputationList <- function( x , subset ,
					select = NULL , expr_subset = NULL , 
					index = NULL , ... ){
		CALL <- match.call()
		if (missing(subset)){  subset <- TRUE	}		
		datlist2 <- subset_datlist( datlist = x , subset = subset ,
					   select = select , expr_subset=expr_subset,
					   index = index , toclass = "imputationList")					   
		datlist2$call <- CALL
		return(datlist2)
				}
#---------------------------------------------------------------	
