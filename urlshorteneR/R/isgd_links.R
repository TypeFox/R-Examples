#' @title Given a long URL, returns a short is.gd link.
#' 
#' @description See \url{https://is.gd/apishorteningreference.php}
#' 
#' @param longUrl - The url parameter is the address that you want to shorten. 
#' 
#' @param shorturl - (optional character) You can specify the shorturl parameter if you'd like to
#' pick a shortened URL instead of having is.gd randomly generate one. These must be between 5 and 30 
#' characters long and can only contain alphanumeric characters and underscores. Shortened URLs 
#' are case sensitive. Bear in mind that a desired short URL might already be taken (this is very 
#' often the case with common words) so if you're using this option be prepared to respond to an 
#' error and get an alternative choice from your app's user.
#' 
#' @param logstats - (optional) Adding the parameter logstats=1 turns on logging of detailed 
#' statistics when the shortened URL you create is accessed. This allows you to see how many 
#' times the link was accessed on a given day, what pages referred people to the link, what 
#' browser visitors were using etc. You can access these stats via the link preview page for 
#' your shortened URL (add a hyphen/dash to the end of the shortened URL to get to it). Creating 
#' links with statistics turned on has twice the "cost" towards our rate limit of other shortened 
#' links, so leave this parameter out of your API call if you don't require statistics on usage. See
#' our usage limits page for more information on this \url{https://is.gd/usagelimits.php}.
#' 
#' @param showRequestURL - show URL which has been build and requested from server. 
#' For debug purposes.
#' 
#' @examples
#' ## asd <- isgd_LinksShorten(longUrl = "http://debil.cz/",showRequestURL = TRUE)
#' 
#' @export
isgd_LinksShorten <- function(longUrl = "", logstats = "0", shorturl = NULL, showRequestURL = FALSE) {
  links_shorten_url <- "http://is.gd/create.php?format=json"
  
  query <- list(url = longUrl, logstats = logstats, shorturl = shorturl)
  
  df_link_shorten <- doRequest("GET", links_shorten_url, queryParameters = query, showURL = showRequestURL)
  
  return(df_link_shorten$shorturl)
  
}
#' @title Expand a short URL to a longer one
#' 
#' @inheritParams isgd_LinksShorten
#' 
#' @description See \url{https://is.gd/apilookupreference.php}
#' 
#' @examples 
#' ### isgd_LinksExpand(shorturl = "http://is.gd/4oIAXJ", showRequestURL = TRUE)
#' 
#' @export
isgd_LinksExpand <- function(shorturl = "", showRequestURL = FALSE) {
  links_expand_url <- "http://is.gd/forward.php?format=json"
  
  query <- list(shorturl = shorturl)
  
  df_link_expand <- doRequest("GET", links_expand_url, queryParameters = query, showURL = showRequestURL)
  
  return(df_link_expand$url)
}

