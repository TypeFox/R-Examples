#' Get Access Token
#'
#' Once you have set your Application ID and Secret via \code{\link{secret_id}}, get the access token.
#' 
#' @return a list with 4 things: access_token, token_type, expires_in, scope
#' @export
#' @references \url{https://developer.clarifai.com/}
#' @examples \dontrun{
#' get_token()
#' }

get_token <- function() {
	
    clarifai_check_auth()

    h <- new_handle()
	handle_setopt(h, customrequest = "POST")
	handle_setform(h, grant_type='client_credentials', client_id=Sys.getenv('ClarifaiId'), client_secret=Sys.getenv('ClarifaiSecret'))
	token_con    <- curl_fetch_memory(url="https://api.clarifai.com/v1/token/", handle=h)
	token_info   <- fromJSON(rawToChar(token_con$content))

	Sys.setenv(ClarifaiToken = token_info$access_token)
	return(invisible(token_info))

}


