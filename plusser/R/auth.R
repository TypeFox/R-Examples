if(getRversion() >= "2.15.1")  utils::globalVariables(c(".gpapikey")) ##to make CRAN happy


##' Sets an API key for the Google+ API
##' 
##' This function sets an API key that is then stored invisibly for
##' \code{plusser} to use when accessing the Google+ API. A warning is issued if
##' URL escaping the api key alters it, as Google should provide you with a
##' HTML-safe API key in the first place.
##' 
##' @param apikey The API key as a character string.
##' @return Returns \code{TRUE} if the key was stored successfully.
##' @export
##' @examples
##' setAPIkey("thisIsInvalid")
setAPIkey <- function(apikey) {
  apikey.enc <- curlEscape(apikey)
  if (apikey.enc != apikey) 
    warning("Your API key has been URL encoded. This should not be neccessary. Check if you used the right key.")
  assign("apikey", apikey, envir=gp)
  return(isTRUE(all.equal(get("apikey", envir=gp),apikey)))
}
