
#############################################
# weighted frequency table
weighted_table <- function( x , w = NULL, props=FALSE ){
	
	#**** vector x
	if ( is.vector(x) ){
	    if ( is.null(w) ){ w <- rep(1,length(x)) }
		x1_u <- sort( unique(x) )
		N1 <- length(x1_u)
		res0 <- rep(NA,N1)
		names(res0) <- x1_u
		for (nn in 1:N1){
			res0[nn] <- sum( ( x == x1_u[nn] ) * w , na.rm=TRUE)
					}	
				}
	#**** matrix x
	if ( is.matrix(x)){
	    if ( is.null(w) ){ w <- rep(1,nrow(x)) }
		x1_u <- sort( unique(x[,1]) )
		N1 <- length(x1_u)
		x2_u <- sort( unique(x[,2]) )
		N2 <- length(x2_u)				
		res0 <- matrix(NA,nrow=N1,ncol=N2)
		rownames(res0) <- x1_u
		colnames(res0) <- x2_u		
		for (nn in 1:N1){
		for (mm in 1:N2){
			res0[nn,mm] <- sum( ( x[,1] == x1_u[nn] ) * ( x[,2] == x2_u[mm] ) * w , na.rm=TRUE)
					}	
				}									
			}				
	#**** calculate proportions?
	if (props){
		res0 <- res0 / sum(res0) 
			}
	
	return(res0)


			}
##############################################			