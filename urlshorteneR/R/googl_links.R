#' @title Expand a short URL to a longer one
#'
#' @seealso See \url{https://developers.google.com/url-shortener/v1/getting_started#shorten}
#' @seealso See \url{https://developers.google.com/url-shortener/v1/url/get}
#'
#' @param shortUrl - The short URL, including the protocol. 
#' @param showRequestURL - show URL which has been build and requested from server. For debug 
#' purposes.
#' @param projection - "FULL" - returns the creation timestamp and all available analytics (default) 
#' OR "ANALYTICS_CLICKS" - returns only click counts OR "ANALYTICS_TOP_STRINGS" - returns only top 
#' string counts (e.g. referrers, countries, etc)
#' 
#' @description For the given short URL, the url.get method returns the corresponding long URL and 
#' the status.
#' 
#' @section Quotas: By default, your registered project gets 1,000,000 requests per day for the URL 
#' Shortener API (\code{https://console.developers.google.com/})
#' 
#' @examples 
#' \dontrun{
#' googl_token <- googl_auth(key = "", secret = "")
#' g1 <- googl_LinksExpand(shortUrl = "http://goo.gl/vM0w4",showRequestURL = TRUE)
#' g4 <- googl_LinksExpand(shortUrl="http://goo.gl/vM0w4",projection = "ANALYTICS_TOP_STRINGS")
#' }
#' 
#' @note Returns a dataframe of expanded short URL and a list of its analytics.
#' 
#' @return id - is the short URL you passed in.
#' @return longUrl - is the long URL to which it expands. Note that longUrl may not be present in 
#' the response, for example, if status is "REMOVED".
#' @return status - is "OK" for most URLs. If Google believes that the URL is fishy, status may be 
#' something else, such as "MALWARE".
#'
#' @export
googl_LinksExpand <- function(shortUrl = "", projection = "FULL", showRequestURL = FALSE) {
  links_expand_url <- "https://www.googleapis.com/urlshortener/v1/url"
  
  query <- list(key = googl_token$credentials$access_token, shortUrl = shortUrl, projection = projection)
  
  # call method from ApiKey.R
  df_link_expand <- doRequest("GET", links_expand_url, "googl", queryParameters = query, showURL = showRequestURL)
  df_link_expand_data_analytics <- df_link_expand$analytics
  df_link_expand$analytics <- NULL
  df_link_expand_data <- list(original_data = data.frame(df_link_expand, stringsAsFactors = FALSE),
                              analytics = df_link_expand_data_analytics)
  
  return(df_link_expand_data)
}


#' @title Given a long URL, returns a short Goo.gl link.
#' 
#' @seealso See \url{https://developers.google.com/url-shortener/v1/url/insert}
#' @seealso See \url{https://developers.google.com/url-shortener/v1/getting_started#shorten}
#' @seealso See \url{http://stackoverflow.com/a/13168073}
#' 
#' @param longUrl - a long URL to be shortened (example: http://betaworks.com/).
#' @param showRequestURL - show URL which has been build and requested from server. For debug 
#' purposes.
#' 
#' @description Given a full URL, returns an goo.gl short URL. The returned resource contains the 
#' short URL and the long URL. Note that the returned long URL may be loosely canonicalized, e.g. 
#' to convert "google.com" into "http://google.com/". See the Authentication 
#' \link{googl_auth} section for more details.
#' 
#' @return id is the short URL that expands to the long URL you provided. If your request includes 
#' an auth token, then this URL will be unique. If not, then it might be reused from a previous 
#' request to shorten the same URL.
#' @return longURL - longUrl is the long URL to which it expands. In most cases, this will be the 
#' same as the URL you provided. In some cases, the server may canonicalize the URL. For instance, 
#' if you pass http://www.google.com, the server will add a trailing slash.
#' 
#' @examples 
#' \dontrun{
#' googl_token <- googl_auth(key = "", secret = "")
#' g2 <- googl_LinksShorten(longUrl = "https://developers.google.com/url-shortener/v1/url/insert")
#' }
#' 
#' @export
googl_LinksShorten <- function(longUrl = "", showRequestURL = FALSE) {
  links_shorten_url <- paste0("https://www.googleapis.com/urlshortener/v1/url?key=", googl_token$credentials$access_token)
  
  resource <- paste0("{'longUrl':", paste0("'",longUrl,"'}"))
  
  # call method from ApiKey.R 
  df_link_shorten <- doRequest("POST", links_shorten_url, "googl", queryParameters = resource, showURL = showRequestURL)
  
  df_link_shorten_data <- data.frame(df_link_shorten, stringsAsFactors = FALSE)
  
  # sapply(df_link_shorten_data, class)
  return(df_link_shorten_data)
}

#' @title Retrieves a list of URLs shortened by the authenticated user.
#'
#' @seealso See \url{https://developers.google.com/url-shortener/v1/url/list}
#' @seealso See \url{https://developers.google.com/url-shortener/v1/getting_started#history}
#'
#' @description returns a paginated list of information about short URLs created by a user,
#' sorted in reverse chronological order. Each returned resource contains the short URL,
#' long URL, creation timestamp, and status.
#'
#' @param showRequestURL - show URL which has been build and requested from server. For debug 
#' purposes.
#' @param projection - an optional (!) information to return : "ANALYTICS_CLICKS" - Returns short
#' URL click counts. OR "FULL" - Returns full analytics information (default)
#'
#' @return totalItems - is an approximate number of items in the user's entire history.
#' @return nextPageToken - is an opaque string you can use to get the next page of history.
#' It looks a lot like an ISO 8601 formatted date right now, but you should not count on that
#' being true. The nextPageToken will be present on all but the last page
#' @return long_url - items contains the list of entries for the first "page" of the user's history,
#' in order of descending creation time. The values for each entry are the same as specified in
#' the Analytics section 
#' \url{https://developers.google.com/url-shortener/v1/getting_started#url_analytics}.
#'
#' @examples
#' \dontrun{
#' googl_token <- googl_auth(key = "", secret = "")
#' googl_UserLinkHistory(showRequestURL = TRUE)
#' googl_UserLinkHistory(projection = "FULL", showRequestURL = TRUE)
#' }
#' 
#' @note Requires that the user authenticates with his google account through OAUTH 2.0 ! Thus no API key is necessary
#'
#' @export
googl_UserLinkHistory <- function(projection = "FULL", showRequestURL = FALSE) {
  
  user_linkHistory_url <- "https://www.googleapis.com/urlshortener/v1/url/history"
  
  query <- list(key = googl_token$credentials$access_token, projection = projection)
  
  df_history <- doRequest("GET", user_linkHistory_url, "googl", queryParameters = query, showURL = showRequestURL)
  
  df_history$items <- data.frame(df_history$items, stringsAsFactors = F)
  
  # sapply(df_history$items, class)
  return(df_history$items)
}
