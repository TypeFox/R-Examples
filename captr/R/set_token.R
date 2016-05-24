#' Sets Application Token
#'
#' Captricity requires an application token to use the API. Get the token from \url{https://shreddr.captricity.com/}.
#' The functions looks for \code{CaptricityToken} in the environment. If it doesn't find it or if change is forced, 
#' it looks for arguments passed in the function. If it fails to find that, it asks for input.
#' 
#' Run this function before anything else.
#'
#' @param app_token Application token. Get these from \url{https://shreddr.captricity.com/developer/}. 
#' @param force Force change the \code{CaptricityToken} stored in the environment
#' 
#' @keywords Sets Application Token
#' @export
#' @references \url{https://shreddr.captricity.com/developer/}
#' @examples \dontrun{
#' set_token("app_token")
#' }

set_token <- 
function(app_token=NULL, force=FALSE) {

    env_id <- Sys.getenv('CaptricityToken')
    
    # If you cannot find CaptricityToken in the environment
    if (identical(env_id, "") | !force) {

    	# First look for arguments passed in the function
	    if (!is.null(app_token)) {
	        Sys.setenv(CaptricityToken = app_token)
	       }

		# Else ask user for the details    
	    else {
    		message("Couldn't find env var CaptricityToken. See ?set_token for more details.")
			message("Please enter your CaptricityToken and press enter:")
		  	pat <- readline(": ")
        	Sys.setenv(CaptricityToken = pat)
	        }
    }

}
