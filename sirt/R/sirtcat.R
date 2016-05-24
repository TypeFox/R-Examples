
######################################################
sirtcat <- function( label , time0 , active ){
	if (active){
		z0 <- time0
		cat( label , "  " )
		z1 <- Sys.time()
		print(z1-z0)
		z0 <- z1 	
		zout <- z0
			} else {
		zout <- NULL
				}
	return(zout)
		}
######################################################	