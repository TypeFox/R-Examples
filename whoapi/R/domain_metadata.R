#'@title Identify a Domain's Search Ranking
#'@description \code{domain_rank} provides users with the ability
#'to identify a domain's "search ranking" - how prominent it is
#'according to common internet benchmarks. Specifically, it provides
#'both Alexa reach and popularity scores, and the Google Pagerank ("pr")
#'ranking, which goes from 0 to 10 (and is represented by -1 if it
#'cannot be found).
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{domain_search}} which specifically
#'looks at the search results a particular domain pulls up
#'in various search engines.
#'
#'@examples
#'token <- whoapi_token("demokey")
#'\dontrun{
#'whoapi_domain <- domain_rank(token, "whoapi.com")
#'}
#'@export
domain_rank <- function(token, domain, ...){
  url <- paste0("&r=ranks&domain=", domain)
  result <- whoapi_query(token, url, ...)
  return(result)
}

#'@title Identify a Domain's Search Results Count
#'@description \code{domain_search} allows you to
#'quickly calculate, in an automated way, how many
#'search results Bing and Google return for a particular
#'domain.
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{domain_rank}} which looks more
#'generally at domain popularity according to Alexa and Google Pagerank
#'scores.
#'
#'@examples
#'token <- whoapi_token("demokey")
#'\dontrun{
#'search_results <- domain_search(token, "whoapi.com")
#'}
#'@export
domain_search <- function(token, domain, ...){
  url <- paste0("&r=searchengines&domain=", domain)
  result <- whoapi_query(token, url, ...)
  return(result)
}

#'@title Retrieve Domain Metadata
#'@description \code{domain_metadata} retrieves
#'information about the content on a domain,
#'specifically the title and metadata description
#'fields from its home page
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@examples
#'token <- whoapi_token("demokey")
#'\dontrun{
#'metadata <- domain_metadata(token, "whoapi.com")
#'}
#'@export
domain_metadata <- function(token, domain, ...){
  url <- paste0("&r=meta&domain=", domain)
  result <- whoapi_query(token, url, ...)
  return(result)
}

#'@title Retrive Domain Location Information
#'@description \code{domain_location} returns geographic
#'information about where a domain - or, specifically,
#'its IP address - is located.
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{whois_info}} for more free-form information,
#'including (potentially) the address of the domain holders.
#'
#'@examples
#'token <- whoapi_token("demokey")
#'\dontrun{
#'location_data <- domain_location(token, "whoapi.com")
#'}
#'@export
domain_location <- function(token, domain, ...){
  url <- paste0("&r=geo&domain=", domain)
  result <- whoapi_query(token, url, ...)
  return(result)
}