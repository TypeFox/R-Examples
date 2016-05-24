mice.impute.2l.groupmean <- function (y, ry, x, type , grmeanwarning=TRUE, ...){  
    if ( ( ncol(x) > 2 ) & grmeanwarning )   warning("\nMore than one variable is requested to be aggregated.\n") 
	clusterx <- paste( x[,type==-2] )
	a1 <- rowsum( x[ , type %in% c(1,2) ] , clusterx ,  na.rm= TRUE )	
	a2 <- rowsum( 1+0*x[ , type %in% c(1,2) ] , clusterx ,  na.rm= TRUE )
    i1 <- match( clusterx  , as.numeric(rownames(a1) ) )
	ximp <- a1[i1,,drop=FALSE]  / a2[i1,,drop=FALSE]
    # calculate aggregated values
	colnames(ximp) <- paste( names(type)[ type %in% c(1,2) ] , names(type)[ type == -2 ] , sep="." )
    return(ximp)
}
