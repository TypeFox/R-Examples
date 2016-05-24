

 
############################################################################
# This function reads fwf files                                            #
read.fwf2 <- function( file , format.full , variables = NULL){
	format <- format.full
    ff <- base::readLines( file )
    ind.ff1 <- c( 1, cumsum(format)[- length(format) ] + 1 )
    ind.ff2 <- cumsum(format)
    I <- length(format)
    n <- length( ff )
    dfr <- data.frame( matrix(0 , nrow= n , ncol=I ) )
    for (ii in 1:I){  
			dfr[,ii ] <- as.numeric( substring( ff , ind.ff1[ii] , ind.ff2[ii] )  ) 
					}
    if (!is.null(variables)){ colnames(dfr) <- variables }
    return(dfr)
    } 
############################################################################