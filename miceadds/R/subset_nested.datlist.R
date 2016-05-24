
#######################################################
subset_nested.datlist <- function( datlist , subset = TRUE ,
				select = NULL , expr_subset = NULL , 
				index_between = NULL , index_within = NULL ,
				toclass = "nested.datlist" , simplify = FALSE ){
		CALL <- match.call()
		
		
		#-----------------------------
		#*** check here for classes
		if ( class(datlist) %in% "NestedImputationList" ){
			datlist <- datlist$imputations		
						}
							
		#*****************************
		# check for expr
		expr <- expr_subset 
		is_expr <- FALSE
		apply_select0 <- FALSE
		pf <- parent.frame()
	    if (!is.null(match.call()$expr)){
		    expr1 <- substitute(expr)
			is_expr <- TRUE
			apply_select0 <- TRUE
					}
																																																																							
		#-----------------------------
		#*** start routine
		NB <- length(datlist)
		NW <- length(datlist[[1]])			
		if ( is.null(index_between) ){
			index_between <- 1:NB
						}		
		if ( is.null(index_within) ){
			index_within <- 1:NW
						}								
						
		IMB <- length(index_between)				
		IMW <- length(index_within)				
				
		if( is.null(select) & ( mean( subset ) == 1 ) ){
			apply_select <- FALSE
						} else {
			apply_select <- TRUE			
						}
		if (apply_select0){ apply_select <- TRUE }
		if ( is.null(select) ){
				select <- colnames(datlist[[1]][[1]])			
						}		
						
							
		# initialize object structure				
		datlist2 <- as.list(1:IMB)
			for (ii in 1:IMB){
				datlist2[[ii]] <- as.list( 1:IMW)
					}
		
		for (ii in 1:IMB){
		for (jj in 1:IMW){
			d1 <- datlist[[ index_between[ii] ]][[ index_within[jj] ]]	
			if (is_expr){
				 # h1 <- with( data=d1 , { expr1 } )
				 subset <- eval(expr1, d1, enclos=pf)
						}
			if (apply_select){
					d1 <- subset( d1 , subset=subset , select=select , drop=FALSE)
							 }
			datlist2[[ii]][[jj]] <- d1
						}
					  }
						
						
		#************
		# create object classes
		#---- class datlist
		if (toclass == "nested.datlist" ){
			datlist2 <- nested.datlist_create(datlist2)
			a2 <- attr(datlist2,"Nimp")
			# simplify within nest
			if ( simplify){
				if ( a2[2] == 1 ){									
					IB <- a2[1]
					datlist3 <- as.list(1:IB)
					for (ii in 1:IB){
						datlist3[[ii]] <- datlist2[[ii]][[1]]
						}
					datlist2 <- datlist_create(datlist3)
					simplify <- FALSE		
							}
					}
			# simplify between nest
			if ( simplify){
				if ( a2[1] == 1){
					datlist2 <- datlist_create( datlist2[[1]] )
							}			
						}
			
			attr(datlist2,"call") <- CALL
								   }	
		#---- class imputationList
		if (toclass == "NestedImputationList" ){
			datlist2 <- NestedImputationList(datlist2)
			datlist2$call <- CALL
			
			a2 <- datlist2$Nimp
			# simplify within nest
			if ( simplify){
				if ( a2[2] == 1 ){									
					IB <- a2[1]
					datlist3 <- as.list(1:IB)
					for (ii in 1:IB){
						datlist3[[ii]] <- datlist2$imputations[[ii]][[1]]
						}
					datlist2 <- mitools::imputationList(datlist3)
					simplify <- FALSE		
							}
					}
			# simplify between nest
			if ( simplify){
				if ( a2[1] == 1){
					datlist2 <- mitools::imputationList( datlist2$imputations[[1]] )
							}			
						}						
								   }

								   
		return(datlist2)
		}
		
		
		
############################################################
# object of class nested.datlist
subset.nested.datlist <- function( x , subset ,
					select = NULL , expr_subset = NULL , 
				index_between = NULL , index_within = NULL , 
				simplify=FALSE , ... ){
		CALL <- match.call()
		if (missing(subset)){  subset <- TRUE	}
		datlist2 <- subset_nested.datlist( datlist = x , subset = subset ,
					   select = select , index_between = index_between , 
					   index_within = index_within , simplify = simplify , 
					   toclass = "nested.datlist")					   
		attr(datlist2,"call") <- CALL
		return(datlist2)
				}
#---------------------------------------------------------------				

#---------------------------------------------------------------						
# object of class imputationList
subset.NestedImputationList <- function( x , subset ,
					select = NULL , expr_subset = NULL , 
				index_between = NULL , index_within = NULL , 
				simplify = FALSE , ... ){
		CALL <- match.call()
		if (missing(subset)){  subset <- TRUE	}		
		datlist2 <- subset_nested.datlist( datlist = x , subset = subset ,
					   select = select , index_between = index_between , 
					   index_within = index_within , simplify = simplify , 
					   toclass = "NestedImputationList")						   
		datlist2$call <- CALL
		return(datlist2)
				}
#---------------------------------------------------------------	
