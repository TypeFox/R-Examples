


print.object.summary <- function( obji , digits ){
    V <- ncol(obji)
    for (vv in 1:V){
    if ( is.numeric( obji[,vv] ) ){ 
		obji[,vv] <- round( obji[,vv] , digits ) 
				}
                }
    print(obji)
        }