
add.colnames.resp <- function(resp){
	if( is.null(colnames(resp)) ){
		I <- ncol(resp)
		colnames(resp) <- paste0("I",1:I) 
					}
	return(resp)
		}