
#' Simple call to test / revalidate token
#' @return list with result == "success"
#' @export
keepalive <- function(){
	
	#test if we can do a read
	oh.user.read("");
	
	#if no error: return success
	return(list(result="success"));	
}
