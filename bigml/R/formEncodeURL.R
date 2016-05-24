#' A simple function to turn named arguments into a form-encoded string
#' @export
#' @param \dots arbitrary named arguments that will become part of a form-encoded url.
#' @param a something
#' @details This function is called in every BigML API function.  It helps
#' 	build the URL that requests are forwarded to.  It automatically adds any
#'	default user and api key settings specified by
#'	\code{\link{setCredentials}}.  However, it also can be used to access
#'	advanced options that are otherwise undocumented here.  For instance, it's
#'	possible to filter and/or sort on a number of different api requests,
#'	using a number of different fields (e.g., see the documentation on
#'	\href{https://bigml.com/developers/datasets#d_list}{listing and sorting
#'	datasets}.)
#'	Other usage includes specifying \code{username} and \code{api_key}
#'	for individual API requests; or \code{limit} or \code{offset} parameters
#' 	useful for paging through list requests.  Finally, it's possible to
#'	enable a simple debug mode by passing debug=TRUE.  This will print the
#'	url request string to the screen, along with any posted json objects.
#' @return form-encoded string result
#' @template author
#' @examples
#' \dontrun{
#' formEncodeURL(username="user1", api_key="test", limit=100, debug=TRUE)
#' # "?username=user1&api_key=test&limit=100&debug=TRUE"
#' }
formEncodeURL <-
function (a,...)
{
    opt_args = list(...)
	sys_user = Sys.getenv("BIGML_USERNAME")
	sys_key = Sys.getenv("BIGML_API_KEY")
	if (!"username" %in% opt_args && sys_user != ''){
		opt_args$username = sys_user
	}
	if (!"api_key" %in% opt_args && sys_key != ''){
		opt_args$api_key = sys_key
	}
    res = ""
    for (v in names(opt_args)) {
        encoded = curlPercentEncode(opt_args[v])
        res = paste(res, v, "=", encoded, "&", sep = "")
    }
	if (nchar(res) <= 0){
		return('')
	} else{
		res = strtrim(res, nchar(res) - 1)
	    return(paste("?", res, sep = ""))
	}

}
