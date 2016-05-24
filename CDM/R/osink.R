
# 	osink( file = file , suffix = "__SUMMARY.Rout" )

#   csink( file = file )


osink <- function( file , suffix){
	if ( ! is.null( file ) ){
		base::sink( paste0( file , suffix) , split=TRUE )
						}
				}
				
csink <- function( file){
	if ( ! is.null( file ) ){  
	   base::sink()	
				}	
					}