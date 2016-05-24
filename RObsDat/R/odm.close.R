odm.close <- function(){
	handler <- options("odm.handler")[[1]]
	if(!is.null(handler)){
		dbDisconnect(handler@con)
		options(odm.handler=NULL)
		cat("Connection closed\n")
	}

}
