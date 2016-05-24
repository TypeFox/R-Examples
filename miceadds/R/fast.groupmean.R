fast.groupmean <- function( data , group , weights=NULL , extend=FALSE){
	groups <- sort( unique( group ) )
	index.group <- match( group , groups )
	if ( is.null(weights) ){
		Ngroup <- rowsum( 1-is.na(data) , index.group ) 
	    data1 <- rowsum( data , index.group , na.rm=TRUE)		
						} else {
		Ngroup <- rowsum( weights*(1-is.na(data)) , index.group )
	    data1 <- rowsum( data*weights , index.group , na.rm=TRUE)		
							}
	colnames(data1) <- colnames(data)
	data1 <- data1 / Ngroup
	data1 <- data.frame( "group" = groups , data1 )
	if (extend){
	   data1 <- data1[ index.group , ]
	   rownames(data1) <- NULL
				}
	return(data1)
            }
#.....
# ARb 2013-10-25
# extend this function to include weights
# and calculation of standard deviation and skewness
#.....
fast.groupsum <- function( data , group , weights=NULL , extend=FALSE){
	groups <- sort( unique( group ) )
	index.group <- match( group , groups )
	if ( is.null(weights) ){
#		Ngroup <- rowsum( 1-is.na(data) , index.group ) 
	    data1 <- rowsum( data , index.group , na.rm=TRUE)		
						} else {
#		Ngroup <- rowsum( weights*(1-is.na(data)) , index.group )
	    data1 <- rowsum( data*weights , index.group , na.rm=TRUE)		
							}
	colnames(data1) <- colnames(data)
#	data1 <- data1 / Ngroup
	data1 <- data.frame( "group" = groups , data1 )
	if (extend){
	   data1 <- data1[ index.group , ]
	   rownames(data1) <- NULL
				}	
	return(data1)
            }