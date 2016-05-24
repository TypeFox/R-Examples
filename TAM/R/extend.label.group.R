


#######################################################
# extend labels arranged in a group
extend.label.group <- function( label.group ){
    str1 <- unlist( strsplit( label.group , split="__" ) )
    M1 <- min( nchar( str1) )
	M2 <- max( nchar( str1) )
    # ii <- 1
	
	# extract types of labels
	str1a <- strsplit( str1[1] , split="" ) 
	str2a <- strsplit( str1[2] , split="" ) 
	
	dfr2 <- matrix( NA , nrow=M2 , ncol=5)
	dfr2[1:M1,1] <- unlist( str1a)
	dfr2[1:M2,2] <- unlist( str2a)
	dfr2 <- as.data.frame(dfr2)
	dfr2[,3] <- paste( dfr2[,1] ) %in% c( LETTERS , letters ) 
	dfr2[,4] <- paste( dfr2[,2] ) %in% c( LETTERS , letters ) 
	dfr2[,5] <- paste(dfr2[,1]) == paste(dfr2[,2])
	ii0 <- 0
	for (ii in 1:M1){
		if ( dfr2[ii,5] ){ 
			ii0 <- ii0+1 
				}	
					}
	if (ii0==M1){ ii0 <- ii0-1 }		
#    for (mm in 1:M1){
#        if ( substring( str1[1] , mm , mm) == substring( str1[2] , mm , mm) ){
#            ii <- ii + 1 
#                } else {
#            break 
#                    }
#            }
	ii <- ii0
    l1 <- as.numeric( substring( str1 , ii+1 )    )
    str2 <- seq( l1[1] , l1[2] )
    str2 <- paste0( substring( str1 , 1 , ii  )[1] , str2 )
	
    return( str2 )
        }
#######################################################

