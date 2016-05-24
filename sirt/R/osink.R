

osink <- function( file , prefix){
	if ( ! is.null( file ) ){
		base::sink( paste0( file , prefix) , split=TRUE )
						}
				}
				
csink <- function( file){
	if ( ! is.null( file ) ){  sink()	}	
					}