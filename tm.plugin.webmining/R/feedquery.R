#' @title Buildup string for feedquery. 
#' @description Function has partly been taken from \code{\link[RCurl]{getForm}} function. 
#' Generally, a feed query is a string built up as follows: \cr
#' \code{<url>?<param1=value1>&<param2=value2>&...&<paramN=valueN>} \cr
#' By specifying a feed url and parameter--value pairs (as list) we can easily
#' generate a feed query in R.
#' @author Mario Annau
#' @param url character specifying feed url
#' @param params list which contains feed parameters, e.g. list(param1="value1", param2="value2")
#' @seealso \code{\link{xmlNode}} \code{\link{getForm}}
#' @examples
#' \dontrun{
#' feedquery(url = "http://dummy.com", 
#' params = list(param1 = "value1", param2 = "value2"))
#' } 
#' @export 
#' @importFrom RCurl curlEscape
feedquery <-
function(url, params){
	els <- lapply(names(params), function(n) {		
		paste(n, curlEscape(params[[n]]), sep = "=")
	})
	names(els) <- names(params)
	
	feeds <- ""
	for(i in names(els)){
		if(feeds[1] == ""){
			sep = ""
		}
		else{
			sep = "&"
		}
		feeds <- paste(feeds, els[[i]], sep = sep)
	}

	feeds <- paste(url, feeds, sep = "?")
	return(feeds)
}