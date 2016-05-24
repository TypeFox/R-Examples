
########################################################
# auxiliary function
ma.wtd.aux.data <- function(data , weights , vars = NULL ){
	is_dfr <- TRUE		# default is class data frame

	
	#---- mids or mids.1chain
	if ( class(data) %in% c("mids","mids.1chain","mids.nmi") ){
		data <- mids2datlist( data )
								}		
	
	#----- NestedImputationList
	if ( class(data) == "NestedImputationList" ){
		is_dfr <- FALSE
		data <- data$imputations
		class(data) <- "nested.datlist"
						}						
	
	
	#--------------------
	# conversion in case of a nested datalist	
	if ( class(data) == "nested.datlist" ){
		is_dfr <- FALSE
		data <- nesteddatlist2datlist(data)
						}											
						
	#---- imputationList
	if ( class(data) == "imputationList" ){
		data <- data$imputations
		class(data) <- "datlist"
								}

	
	#--------------------
	# conversion in case of a datlist
	if ( class(data) == "datlist" ){
		is_dfr <- FALSE
		if ( ! is.null(vars) ){
			M <- length(data)
			for (ii in 1:M){
			    dat0 <- data[[ii]]
				data[[ii]] <- dat0[ , vars, drop=FALSE ]		
								}
							}
						}
					
						
	#-------------------					
	# conversion in case of a data frame
	if ( is_dfr ){
	    data0 <- data
		if ( ! is.null(vars) ){
			data0 <- data0[ , vars , drop=FALSE ]
								}			
		data <- list(1)
		data[[1]] <- data0
				}
	#-------------------
	# creation of weights (if needed)
	if ( is.null(weights) ){
		weights <- rep(1 , nrow(data[[1]] ) )
						}
	res <- list( data = data , weights = weights )
			}
#########################################################	



########################################################
# weighted mean
ma.wtd.meanNA <- function( data , weights = NULL  , vars = NULL ){    
	#*** pre-processing
	res <- ma.wtd.aux.data(data=data , weights=weights , vars = vars )
	data <- res$data
	weights <- res$weights	
	M <- length(data)
    #*** weighted means
	res <- matrix( NA , nrow=M , ncol= ncol(data[[1]]) )	
	for (ii in 1:M){
		data1 <- data[[ii]]
		dataResp <- 1 - is.na( data1 )
		data1[ is.na(data1) ] <- 0
		data1 <- as.matrix( data1 )        
		# calculate means
		sumweight <- colSums( dataResp * weights )
		M_vars <- colSums( data1 *  weights ) / sumweight
		res[ii,] <- M_vars
					}	
	res <- colMeans(res)
	names(res) <- colnames(data[[1]])
    return( res )
        }
#########################################################



###############################################################################
# weighted SD
ma.wtd.sdNA <- function( data , weights = NULL , vars = NULL , method="unbiased" ){
	#*** pre-processing
	res <- ma.wtd.aux.data(data=data , weights=weights , vars = vars )
	data <- res$data
	weights <- res$weights	
	M <- length(data)
    #*** weighted sds
	res <- matrix( NA , nrow=M , ncol= ncol(data[[1]]) )	
	for (ii in 1:M){		
        data1 <- data[[ii]]	
		dataResp <- 1 - is.na( data1 )
		data1[ is.na(data1) ] <- 0
		data1 <- as.matrix( data1 )        
		# calculate means
		sumweight <- colSums( dataResp * weights )
		M_vars <- colSums( data1 *  weights ) / sumweight
		M_varsM <- matrix( M_vars , nrow= nrow(data1) , ncol=length(M_vars ) , byrow=TRUE )
		data1adj <- ( data1 - M_varsM ) * dataResp # take care of missings    
		w1 <- colSums( dataResp * weights )
		sdx <- sqrt( colSums( data1adj^2 * weights ) /  w1 )
		# adjustment of covariance
		if (method == "unbiased" ){
	        wgtadj <- w1 - colSums( dataResp * weights^2 ) / w1		
			wgtadj <- w1 / wgtadj
			# wtd.var ... 
			#     sum(weights * ((x - xbar)^2))/(sw - sum(weights^2)/sw)	
			sdx <- sqrt(wgtadj) * sdx	
								}		
		res[ii,] <- sdx
					}
	res <- colMeans(res)
	names(res) <- colnames(data[[1]])
    return( res )
        }
###############################################################################




###############################################################################
# weighted covariance
ma.wtd.covNA <- function( data , weights = NULL , vars = NULL ,  method="unbiased" ){
	#*** pre-processing
	res <- ma.wtd.aux.data(data=data , weights=weights , vars = vars )
	data <- res$data
	weights <- res$weights	
	M <- length(data)
    #*** weighted covariance
	V <- ncol(data[[1]])	
	res <- array( NA , dim=c(M,V,V) )
	dimnames(res)[[2]] <- colnames(data[[1]] )
	dimnames(res)[[3]] <- colnames(data[[1]] )
	for (ii in 1:M){
		data1 <- data[[ii]]
		if ( ! is.null(vars) ){
			data1 <- data1[ , vars , drop=FALSE ]
								}	
		dataResp <- 1 - is.na( data1 )
		data1[ is.na(data1) ] <- 0
		data1 <- as.matrix( data1 )        
		# calculate means
		sumweight <- colSums( dataResp * weights )
		M_vars <- colSums( data1 *  weights ) / sumweight
		M_varsM <- matrix( M_vars , nrow= nrow(data1) , ncol=length(M_vars ) , byrow=TRUE )
		data1adj <- ( data1 - M_varsM ) * dataResp # take care of missings    
		sqrtweights <- sqrt( weights )
		# calculate weighted covariance
		# cross-products
		covXY <- crossprod( data1adj * sqrtweights )
		w1 <- covWXY <- crossprod( dataResp * sqrtweights )
		covXY <- covXY / covWXY
		# adjustment of covariance
		if (method == "unbiased" ){
	        # wgtadj <- w1 - colSums( dataResp * weights^2 ) / w1		
			# wgtadj <- w1 / wgtadj		
			wgtadj <- crossprod( dataResp * weights^2 )
			wgtadj <- ( covWXY^2 - crossprod( dataResp * weights^2 ) ) / covWXY^2
			wgtadj <- 1 / wgtadj
			covXY <- wgtadj * covXY
								}						
		res[ii,,] <- covXY
					}
	    res <- colMeans(res)
		return( res )
        }
###############################################################################
# weighted correlation using cov2cor function
ma.wtd.corNA <- function( data , weights = NULL , vars = NULL , method="unbiased" ){
	#*** pre-processing
	res <- ma.wtd.aux.data(data=data , weights=weights , vars = vars )
	data <- res$data
	weights <- res$weights	
	M <- length(data)
    #*** weighted covariance
	V <- ncol(data[[1]])	
	res <- array( NA , dim=c(M,V,V) )
	dimnames(res)[[2]] <- colnames(data[[1]] )
	dimnames(res)[[3]] <- colnames(data[[1]] )	
	for (ii in 1:M){
		data1 <- data[[ii]]
		if ( ! is.null(vars) ){
			data1 <- data1[ , vars , drop=FALSE ]
								}	
		dataResp <- 1 - is.na( data1 )
		data1[ is.na(data1) ] <- 0
		data1 <- as.matrix( data1 )        
		# calculate means
		sumweight <- colSums( dataResp * weights )
		M_vars <- colSums( data1 *  weights ) / sumweight
		M_varsM <- matrix( M_vars , nrow= nrow(data1) , ncol=length(M_vars ) , byrow=TRUE )
		data1adj <- ( data1 - M_varsM ) * dataResp # take care of missings    
		sqrtweights <- sqrt( weights )
		# calculate weighted covariance
		# cross-products
		covXY <- crossprod( data1adj * sqrtweights )
		covWXY <- crossprod( dataResp * sqrtweights )
		covXY <- covXY / covWXY
		# adjustment of covariance
		if (method == "unbiased" ){
			wgtadj <- crossprod( dataResp * weights^2 )
			wgtadj <- ( covWXY^2 - crossprod( dataResp * weights^2 ) ) / covWXY^2
			wgtadj <- 1 / wgtadj
			covXY <- wgtadj * covXY
								}					
		res[ii,,] <- stats::cov2cor( covXY )
					}
	    res <- colMeans(res)
		return( res )
			}
################################################################			