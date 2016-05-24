#' @title Query for a Bitlink based on a long URL.
#' 
#' @description See \url{http://dev.bitly.com/links.html#v3_link_lookup}
#'
#' @param url - one long URLs to lookup.
#' @param showRequestURL - show URL which has been build and requested from server. For debug 
#' purposes.
#'      
#' @return url - an echo back of the url parameter.
#' @return aggregate_link - the corresponding bitly aggregate link (global hash).
#' 
#' @examples
#' \dontrun{ 
#' bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8", 
#' secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")
#' bitly_LinksLookup(url = "http://www.seznam.cz/")
#' bitly_LinksLookup(url = "http://www.seznam.cz/", showRequestURL = TRUE) 
#'
#' manyUrls <- list("http://www.seznam.cz/", "http://www.seznamasdas.cz/", 
#' "http://www.seznam.cz/asadasd", "http://www.seznam.cz/adqwrewtregt")
#' for (u in 1:length(manyUrls)) {
#'    print(bitly_LinksLookup(url = manyUrls[[u]], showRequestURL = TRUE))
#' }
#' }
#' 
#' @export
bitly_LinksLookup <- function(url, showRequestURL = FALSE) {
  
  links_lookup_url <- "https://api-ssl.bitly.com/v3/link/lookup"
  
  query <- list(access_token = bitly_token$credentials$access_token, url = url)
  
  # call method from ApbiKey.R
  df_link_lookup <- doRequest("GET", links_lookup_url, "bitly", query, showURL = showRequestURL)
  df_link_lookup_data <- df_link_lookup$data$link_lookup
  
  # sapply(df_link_lookup_data, class)
  return(df_link_lookup_data)
}

#' @title Used to return the page title for a given Bitlink.
#' 
#' @description See \url{http://dev.bitly.com/links.html#v3_info}
#' 
#' @note Either shortUrl or hash must be given as a parameter (or both).
#' @note The maximum number of shortUrl and hash parameters is 15.
#' 
#' @param hashIN - refers to one bitly hashes, (e.g.:  2bYgqR or a-custom-name). Required
#' @param shortUrl - refers to one Bitlinks e.g.: http://bit.ly/1RmnUT or http://j.mp/1RmnUT. 
#' Optional.
#' @param expand_user - optional true|false (default) - include extra user info in response.
#' @param showRequestURL - show URL which has been build and requested from server. 
#' For debug purposes.
#'  
#' @return short_url - this is an echo back of the shortUrl request parameter.
#' @return hash - this is an echo back of the hash request parameter.
#' @return user_hash - the corresponding bitly user identifier.
#' @return global_hash - the corresponding bitly aggregate identifier.
#' @return error - indicates there was an error retrieving data for a given shortUrl or hash. An 
#' example error is "NOT_FOUND".
#' @return title - the HTML page title for the destination page (when available).
#' @return created_by - the bitly username that originally shortened this link, if the 
#' link is public. Otherwise, null.
#' @return created_at - the epoch timestamp when this Bitlink was created.
#' 
#' @examples
#' \dontrun{ 
#' bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8", 
#' secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")
#' bitly_LinksInfo(shortUrl = "http://bit.ly/DPetrov")
#' bitly_LinksInfo(hash = "DPetrov", showRequestURL = TRUE) 
#' bitly_LinksInfo(hash = "DPetrov", expand_user = "true")
#' 
#' ## hash is the one which is only returned. Dont use
#' bitly_LinksInfo(shortUrl = "on.natgeo.com/1bEVhwE", hash = "DPetrov") 
#' 
#' ## manyHashes <- list("DPetrov", "1QU8CFm", "1R1LPSE", "1LNqqva")
#' ## for (u in 1:length(manyHashes)) {
#' ##   print(bitly_LinksInfo(hashIN = manyHashes[[u]], showRequestURL = TRUE))
#' ## }
#' }
#' 
#' @export
bitly_LinksInfo <- function(hashIN = NULL, shortUrl = NULL, expand_user = "true", showRequestURL = FALSE) {
  
  links_info_url <- "https://api-ssl.bitly.com/v3/info"
  
  if (is.null(hashIN)) {
    query <- list(access_token = bitly_token$credentials$access_token, shortUrl = shortUrl, expand_user = expand_user)
  } else {
    query <- list(access_token = bitly_token$credentials$access_token, hash = hashIN, expand_user = expand_user)
  }
  
  # call method from ApbiKey.R
  df_link_info <- doRequest("GET", links_info_url, "bitly", query, showURL = showRequestURL)
  
  df_user_info_data <- data.frame(t(sapply(unlist(df_link_info$data$info), c)), stringsAsFactors = FALSE)
  df_user_info_data$created_at <- as.POSIXct(as.integer(df_user_info_data$created_at), 
                                             origin = "1970-01-01", tz = "UTC")
  
  # sapply(df_link_info_data, class)
  return(df_user_info_data)
}


#' @title Given a bitly URL or hash (or multiple), returns the target (long) URL.
#' 
#' @description See \url{http://dev.bitly.com/links.html#v3_expand}
#'
#' @param hashIN - refers to one bitly hashes, (e.g.:  2bYgqR or a-custom-name). Required
#' @param shortUrl - refers to one Bitlinks e.g.: http://bit.ly/1RmnUT or http://j.mp/1RmnUT. 
#' Optional.
#' @param showRequestURL - show URL which has been build and requested from server. 
#' For debug purposes.
#' 
#' @section TODO: or more URLs  Up TO 15
#' 
#' @note Either shortUrl or hash must be given as a parameter.
#' @note The maximum number of shortUrl and hash parameters is 15.
#' 
#' @return short_url - this is an echo back of the shortUrl request parameter.
#' @return hash - this is an echo back of the hash request parameter.
#' @return user_hash - the corresponding bitly user identifier.
#' @return global_hash - the corresponding bitly aggregate identifier.
#' @return error - indicates there was an error retrieving data for a given shortUrl or hash. An 
#' example error is "NOT_FOUND".
#' @return long_url - the URL that the requested short_url or hash points to.
#' 
#' @examples
#' \dontrun{ 
#' bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8", 
#' secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")
#' bitly_LinksExpand(shortUrl = "http://bit.ly/DPetrov")
#' bitly_LinksExpand(hash = "DPetrov", showRequestURL = TRUE) 
#' bitly_LinksExpand(hash = "DPetrov")
#' bitly_LinksExpand(shortUrl = "on.natgeo.com/1bEVhwE", hash = "1bEVhwE")
#' 
#' ## manyHashes <- list("DPetrov", "1QU8CFm", "1R1LPSE", "1LNqqva")
#' ## for (u in 1:length(manyHashes)) {
#' ##   print(bitly_LinksExpand(hashIN = manyHashes[[u]], showRequestURL = TRUE))
#' ## }
#' }
#' 
#' @export
bitly_LinksExpand <- function(hashIN = NULL, shortUrl = NULL, showRequestURL = FALSE) {
  
  links_expand_url <- "https://api-ssl.bitly.com/v3/expand"
  
  if (is.null(hashIN)) {
    query <- list(access_token = bitly_token$credentials$access_token, shortUrl = shortUrl)
  } else {
    query <- list(access_token = bitly_token$credentials$access_token, hash = hashIN)
  }
  
  # call method from ApbiKey.R
  df_link_expand <- doRequest("GET", links_expand_url, "bitly", query, showURL = showRequestURL)
  
  df_link_expand_data <- data.frame(t(sapply(unlist(df_link_expand$data$expand), c)), stringsAsFactors = FALSE)
  
  # sapply(df_link_expand_data, class)
  return(df_link_expand_data)
}


#' @title Given a long URL, returns a short Bit.ly link.
#'
#' @description See \url{http://dev.bitly.com/rate_limiting.html} and 
#' \url{http://dev.bitly.com/links.html#v3_shorten}
#'
#' @param longUrl - a long URL to be shortened (example: http://betaworks.com/).
#' @param domain - (optional) the short domain to use; either bit.ly, j.mp, or bitly.com or 
#' a custom short domain. The default for this parameter is the short domain selected by each 
#' user in their bitly account settings. Passing a specific domain via this parameter will override
#' the default settings.
#' @param showRequestURL - show URL which has been build and requested from server. For debug 
#' purposes.
#' 
#' @note Look in the vignette for bulk shortening of URLs. Each call of this function == 1 API call. 
#' Take that into consideration due to limits etc. 
#' @note The bitly API does not support shortening more than one long URL with a single API call. 
#' Meaning 1 Long URL = 1 Function call.
#' @note Long URLs should be URL-encoded. You can not include a longUrl in the request 
#' that has &, ?, #, or other reserved parameters without first encoding it.
#' @note The default value for the domain parameter is selected by each user from within their bitly 
#' account settings at \url{https://bitly.com/a/settings/advanced}.
#' @note Long URLs should not contain spaces: any longUrl with spaces will be rejected. All spaces 
#' should be either percent encoded %20 or plus encoded +. Note that tabs, newlines and trailing 
#' spaces are all indications of errors. Please remember to strip leading and trailing whitespace 
#' from any user input before shortening.
#' 
#' @return new_hash - designates if this is the first time this long_url was shortened by this user. 
#' The return value will equal 1 the first time a long_url is shortened. It will also then be added 
#' to the user history.
#' @return hash - a bitly identifier for long_url which is unique to the given account.
#' @return long_url - an echo back of the longUrl request parameter. This may not always be equal to 
#' the URL requested, as some URL normalization may occur (e.g., due to encoding differences, or case 
#' differences in the domain). This long_url will always be functionally identical the the request 
#' parameter. 
#' @return global_hash - a bitly identifier for long_url which can be used to track aggregate stats 
#' across all Bitlinks that point to the same long_url.
#' @return url - the actual Bitlink that should be used, and is a unique value for the given Bitly 
#' account.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8", 
#' secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")
#' bitly_LinksShorten(longUrl = "http://slovnik.seznam.cz/")
#' bitly_LinksShorten(longUrl = "https://travis-ci.org/dmpe/rbitly/builds/68231423",domain = "j.mp")
#' }
#' 
#' 
#' @export
bitly_LinksShorten <- function(longUrl, domain = NULL, showRequestURL = FALSE) {
  
  links_shorten_url <- "https://api-ssl.bitly.com/v3/shorten"
  
  query <- list(access_token = bitly_token$credentials$access_token, longUrl = longUrl, domain = domain)
  
  # call method from ApiKey.R
  df_link_shorten <- doRequest("GET", links_shorten_url, "bitly", query, showURL = showRequestURL)
  
  df_link_shorten_data <- data.frame(t(sapply(df_link_shorten$data, c)), stringsAsFactors = FALSE)
  
  # sapply(df_link_shorten_data, class)
  return(df_link_shorten_data)
}
