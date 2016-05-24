#'@title Identify if a Domain is Taken
#'@description \code{is_taken} identifies if a domain is taken - if
#'it's been reserve or registered.
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@return a boolean TRUE or FALSE indicating, respectively, whether the domain
#'is reserved/registered, or not.
#'
#'@seealso \code{\link{whois_info}} for more information about a domain's status,
#'including when the registration expires and who has registered it.
#'
#'@examples
#'#Check if whoapi.com is taken
#'token <- whoapi_token("demokey")
#'\dontrun{
#'is_taken(token, "whoapi.com")
#'}
#'#[1] TRUE
#'
#'@export
is_taken <- function(token, domain, ...){
  url <- paste0("&r=taken&domain=", domain)
  result <- whoapi_query(token, url, ...)$taken
  return(char_convert(result))
}

#'@title Identify if a Domain is Blacklisted
#'@description \code{is_blacklisted} checks whether a domain is on prominent
#'spam blacklists (or not).
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@return a list containing a boolean value, "blacklisted", which identifies whether a
#'domain is blacklisted by \emph{any} of the checked services, followed by a breakdown of what the
#'status is on each particular service.
#'
#'@export
is_blacklisted <- function(token, domain, ...){
  url <- paste0("&r=blacklist&domain=", domain)
  result <- whoapi_query(token, url, ...)
  result$status <- NULL
  result$blacklisted <- char_convert(result$blacklisted)
  result$blacklists <- lapply(result$blacklists, function(x){
    x$blacklisted <- char_convert(x$blacklisted)
    return(x)
  })
  return(result)
}

#'@title Get Registration Information for a Domain
#'@description \code{whois_info} grabs detailed information about a domain's registration,
#'including (but not limited to) who it's registered to, what its status is, when it
#'was registered and when it expires.
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{is_taken}} for simpy determining if a domain is registered
#'at all, or \code{\link{certificate_info}} for information about a domain's
#'SSL certificates.
#'
#'@examples
#'token <- whoapi_token("demokey")
#'\dontrun{
#'whoapi_whois_info <- whois_info(token, "whoapi.com")
#'}
#'@export
whois_info <- function(token, domain, ...){
  url <- paste0("&r=whois&domain=", domain)
  result <- whoapi_query(token, url, ...)
  return(result)
}

#'@title Get Certificate Information for a Domain
#'@description \code{certificate_info} retrieves information
#'about a domain's SSL certificate (if present).
#'
#'@param token a token generated with \code{\link{whoapi_token}}
#'
#'@param domain a domain name
#'
#'@param ... further arguments to pass to httr's GET.
#'
#'@seealso \code{\link{whois_info}} for information about the domain
#'in general.
#'
#'@examples
#'token <- whoapi_token("demokey")
#'\dontrun{
#'whoapi_cert <- certificate_info(token, "whoapi.com")
#'}
#'@export
certificate_info <- function(token, domain, ...){
  url <- paste0("&r=cert&domain=", domain)
  result <- whoapi_query(token, url, ...)
  return(result)
}