#' Get Information
#'
#' Basic information about the application --- what are maximum sizes allowed, 
#' default language, max. and min. image and video size, max. batch size, etc.
#' 
#' @return Named list of length 3: \code{status_code}, \code{status_msg}, and \code{results}. 
#' \code{results} is a named list of length 12. Contains information about max. and 
#' min. image and video size allowed etc. 
#' 
#' Prints \code{status_msg} by default 
#' 
#' @export
#' @references \url{https://developer.clarifai.com/}
#' @examples \dontrun{
#' get_info()
#' }

get_info <- function() {

	clarifai_check_token()
		
    h <- new_handle()
	handle_setopt(h,  customrequest = "GET")
	handle_setheaders(h, "Authorization" = paste0("Bearer ", Sys.getenv("ClarifaiToken")))
	info_con   <- curl_fetch_memory("https://api.clarifai.com/v1/info/", handle=h)
	info       <- fromJSON(rawToChar(info_con$content))

	# Print some important things
	cat("Status message: ", info$status_msg, "\n")

	return(invisible(info))

}