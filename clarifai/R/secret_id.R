#' Sets Application ID and Secret
#'
#' Set Client ID and Secret. Needed for interfacing with Clarifai. Run this before anything else.
#' 
#' The function looks for \code{ClarifaiId} and \code{ClarifaiSecret} in the environment. If it doesn't find them and if we don't want to force
#' change in them, it looks for arguments. And if no arguments are passed, it asks for user to input the values. 
#' 
#' @param appdetails A vector of client_id, client_secret. Get these from \url{https://developer.clarifai.com/}. 
#' Set them before you use other functions.
#' @param force force reset client id and secret
#' @keywords Sets Client ID and Secret
#' @export
#' @references \url{https://developer.clarifai.com/}
#' @examples \dontrun{
#' setapp(c("client_id", "client_secret"))
#' }

secret_id <- function(appdetails=NULL, force=FALSE) {

    env_id <- Sys.getenv('ClarifaiId')
    env_pass <- Sys.getenv('ClarifaiSecret')
    
    # If you cannot find ClarifaiId or ClarifaiSecret in the environment
    if ((identical(env_id, "") | identical(env_pass, "")) | !force) {

    	# First look for arguments passed in the function
	    if (!is.null(appdetails)) {
	        Sys.setenv(ClarifaiId = appdetails[1])
	        Sys.setenv(ClarifaiSecret = appdetails[2])
	       }

		# Else ask user for the details    
	    else {
    		message("Couldn't find env var ClarifaiId or ClarifaiSecret. See ?setapp for more details.")
			message("Please enter your ClarifaiId and press enter:")
		  	pat <- readline(": ")
        	Sys.setenv(ClarifaiId = pat)
        	message("Now please enter your ClarifaiSecret and press enter:")
		  	pat <- readline(": ")
        	Sys.setenv(ClarifaiSecret = pat)
	        }
    }
}

   