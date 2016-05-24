#' @title Returns numbers of click on a link.
#' 
#' @description Returns the number of clicks on a single Bitlink of the authenticated user.
#' 
#' @seealso \url{http://dev.bitly.com/link_metrics.html#v3_link_clicks}
#' 
#' @param link - a Bitlink.
#' @param limit - 1 to 1000 (default=1000).
#' @param rollup - true (default) or false.  Return data for multiple units rolled up to a 
#' single result instead of a separate value for each period of time.
#' @param units - an integer representing the time units to query data for. Pass -1 to return all 
#' units of time.
#' @param unit - minute, hour, day, week or month, default: day; Note: when unit is minute the 
#' maximum value for units is 60. 
#' @param showRequestURL - show URL which has been build and requested from server. For debug 
#' purposes.
#' 
#' @return clicks - the number of clicks on the specified Bitlink.
#' @return dt - time in UTC format (only when rollup = "false")
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsClicks(link = "http://bit.ly/DPetrov", unit = "day", units = -1, limit = 100)
#' }
#' 
#' @export
bitly_LinksMetricsClicks <- function(link, limit = 1000, unit = c("minute", "hour", "day", "week", "month"),
                                units = -1, rollup = "true", showRequestURL = FALSE) {
  unit_matched <- match.arg(unit)
  
  link_metrics_clicks_url <- "https://api-ssl.bitly.com/v3/link/clicks"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link, limit = limit, unit = unit_matched, units = units, 
                rollup = rollup)
  
  # call method from ApiKey.R
  df.link_metrics_clicks <- doRequest("GET", link_metrics_clicks_url, service = "bitly", query, showURL = showRequestURL)
  
  df.link_metrics_clicks_data <- df.link_metrics_clicks$data$link_clicks
  
  if (rollup == "true") {
    return(df.link_metrics_clicks_data)
    
  } else {
    df.link_metrics_clicks_data$dt <- as.POSIXct(df.link_metrics_clicks_data$dt, origin = "1970-01-01", tz = "UTC")
  }
  
  # sapply(df.link_metrics_clicks_data, class)
  
  return(df.link_metrics_clicks_data)
  
}


#' @title Returns metrics about the countries from a link.
#' 
#' @description Returns metrics about the countries referring click traffic to a single Bitlink.
#' 
#' @seealso \url{http://dev.bitly.com/link_metrics.html#v3_link_countries}
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @return clicks - the number of clicks referred from this country.
#' @return country - the two-letter code of the referring country.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsCountries(link = "http://bit.ly/DPetrov", unit = "day", units = -1, limit = 100)
#' bitly_LinksMetricsCountries(link = "http://bit.ly/DPetrov", unit = "day", units = -1, limit = 100, 
#' showRequestURL= TRUE)
#' }
#' 
#' @export
bitly_LinksMetricsCountries <- function(link, limit = 1000, unit = c("minute", "hour", "day", "week", "month"),
                                   units = -1, showRequestURL = FALSE) {
  unit_matched <- match.arg(unit)
  
  link_metrics_countries_url <- "https://api-ssl.bitly.com/v3/link/countries"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link, service = "bitly", 
                limit = limit, unit = unit_matched, units = units)
  
  # call method from ApiKey.R
  df_link_metrics_countries <- doRequest("GET", link_metrics_countries_url, "bitly", query, showURL = showRequestURL)
  
  df_link_metrics_countries_data <- df_link_metrics_countries$data$countries
  
  # sapply(df_link_metrics_countries_data, class)
  return(df_link_metrics_countries_data)
}

#' @title Returns users who have encoded this long URL.
#' 
#' @description Returns users who have encoded this long URL (optionally only those in 
#' the requesting user's social graph).
#' 
#' @seealso \url{http://dev.bitly.com/link_metrics.html#v3_link_encoders}
#' 
#' @note Some users may not be returned from this call depending on Bitlink privacy settings.
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @param my_network - true or false (default) restrict to my network.
#' @param subaccounts - (only available to enterprise accounts) false (always default) restrict to 
#' this enterprise account and its subaccounts
#' @param expand_user - true or false (default) include display names of encoders.
#' 
#' @return entries - a mapping of link, user, and ts (when the Bitlink was created) and possible 
#' more depending on input parameters.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsEncoders(link = "http://bit.ly/DPetrov")
#' bitly_LinksMetricsEncoders("http://bit.ly/DPetrov", expand_user = "true", my_network = "false")
#' }
#' 
#' @export
bitly_LinksMetricsEncoders <- function(link, my_network = "false", limit = 25, expand_user = "false", 
                                  subaccounts = "false", showRequestURL = FALSE) {
  
  link_metrics_encoders_url <- "https://api-ssl.bitly.com/v3/link/encoders"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link,  
                limit = limit, my_network = my_network, 
                expand_user = expand_user, subaccounts = subaccounts)
  
  # call method from ApiKey.R
  df_link_metrics_encoders <- doRequest("GET", link_metrics_encoders_url, query, service = "bitly", showURL = showRequestURL)
  
  df_link_metrics_encoders_data <- df_link_metrics_encoders$data$entries
  
  df_link_metrics_encoders_data$ts <- as.POSIXct(as.integer(df_link_metrics_encoders_data$ts),
                                                 origin = "1970-01-01", tz = "UTC")
  # sapply(df_link_metrics_encoders_data, class)
  return(df_link_metrics_encoders_data)
}


#' @title Returns the number of users who have shortened a link.
#' 
#' @description Returns the number of users who have shortened (encoded) a single Bitlink.
#' 
#' @seealso \url{http://dev.bitly.com/link_metrics.html#v3_link_encoders_count}
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @return aggregate_link - the aggregate (global) Bitlink for the provided Bitlink.
#' @return count - the number of bitly users who have shortened (encoded) this link.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsEncodersCount(link = "http://bit.ly/DPetrov")
#' }
#' 
#' @export
bitly_LinksMetricsEncodersCount <- function(link, showRequestURL = FALSE) {
  
  link_metrics_encoders_count_url <- "https://api-ssl.bitly.com/v3/link/encoders_count"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link)
  
  # call method from ApiKey.R
  df_link_metrics_encoders_count <- doRequest("GET", link_metrics_encoders_count_url, service = "bitly", 
                                              query, showURL = showRequestURL)
  df_link_metrics_encoders_count_data <- df_link_metrics_encoders_count$data
  
  # https://stackoverflow.com/questions/4227223/r-list-to-data-frame
  df_link_metrics_encoders_count_data <- data.frame(t(sapply(df_link_metrics_encoders_count_data, c)),  
                                                    stringsAsFactors = FALSE)
  df_link_metrics_encoders_count_data$count <- as.numeric(df_link_metrics_encoders_count_data$count)
  
  # sapply(df_link_metrics_encoders_count_data, class)
  return(df_link_metrics_encoders_count_data)
}

#' @title Returns users who have encoded this link.
#' 
#' @description Returns users who have encoded this link (optionally only those in the requesting 
#' user's social graph), sorted by the number of clicks on each encoding user's link.
#' 
#' @seealso See \url{http://dev.bitly.com/link_metrics.html#v3_link_encoders_by_count}
#' 
#' @note The response will only contain users whose links have gotten at least one click, and 
#' will not contain any users whose links are private.
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @param my_network - true or false (default) restrict to my network
#' @param subaccounts - (only available to enterprise accounts) false (always default) restrict to 
#' this enterprise account and its subaccounts
#' @param expand_user - false (always default) include display names of encoders
#'
#' @return entries - a mapping of link, user, and ts (when the Bitlink was created) and possible 
#' more depending on input parameters.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsEncodersByCount("http://bit.ly/DPetrov", my_network = "false", limit = 100) 
#' bitly_LinksMetricsEncodersByCount("http://bit.ly/DPetrov", my_network = "false", limit = 100, 
#' expand_user = "true")
#' }
#' 
#' @export
bitly_LinksMetricsEncodersByCount <- function(link, limit = 100, my_network = "false", expand_user = "false", 
                                           subaccounts = "false", showRequestURL = FALSE) {
  
  link_metrics_encoders_by_count_url <- "https://api-ssl.bitly.com/v3/link/encoders_by_count"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link, limit = limit, my_network = my_network, 
                expand_user = expand_user, subaccounts = subaccounts)
  
  # call method from ApiKey.R
  df_link_metrics_encoders_by_count <- doRequest("GET", link_metrics_encoders_by_count_url, service = "bitly", 
                                                 query, showURL = showRequestURL)
  
  df_link_metrics_encoders_by_count_data <- data.frame(df_link_metrics_encoders_by_count$data$encoders_by_count)
  df_link_metrics_encoders_by_count_data$ts <- as.POSIXct(as.integer(df_link_metrics_encoders_by_count_data$ts), 
                                                          origin = "1970-01-01", tz = "UTC")
  
  # sapply(df_link_metrics_encoders_by_count_data, class)
  return(df_link_metrics_encoders_by_count_data)
}


#' @title Returns metrics about the domains referring click traffic to a link.
#' 
#' @description Returns metrics about the domains referring click traffic to a single Bitlink.
#' 
#' @seealso \url{http://dev.bitly.com/link_metrics.html#v3_link_referring_domains}
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @return clicks - the number of clicks referred from this domain.
#' @return domain - the domain referring clicks.
#' @return url - the complete URL of the domain referring clicks.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsReferringDomains("http://bit.ly/DPetrov", unit = "day", units=-1, limit = 100)
#' }
#' 
#' @export
bitly_LinksMetricsReferringDomains <- function(link, limit = 1000, unit = c("minute", "hour", "day", "week", "month"), 
                                           units = -1, showRequestURL = FALSE) {
  unit_matched <- match.arg(unit)
  
  link_metrics_referring_domains_url <- "https://api-ssl.bitly.com/v3/link/referring_domains"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link,  
                limit = limit, unit = unit_matched, units = units)
  
  # call method from ApiKey.R
  df_link_metrics_referring_domains <- doRequest("GET", link_metrics_referring_domains_url, service = "bitly", 
                                                 query, showURL = showRequestURL)
  df_link_metrics_referring_domains_data <- df_link_metrics_referring_domains$data$referring_domains
  
  # sapply(df_link_metrics_referring_domains_data, class)
  
  return(df_link_metrics_referring_domains_data)
}


#' @title Returns metrics about the pages referring click traffic to a single Bitlink.
#' 
#' @description \url{http://dev.bitly.com/link_metrics.html#v3_link_referrers}
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @return clicks - the number of clicks referred from this domain.
#' @return referrer - the URL referring clicks.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsReferrers(link = "http://bit.ly/DPetrov",unit = "day", units = -1, limit = 100)
#' }
#' 
#' @export
bitly_LinksMetricsReferrers <- function(link, limit = 1000, unit = c("minute", "hour", "day", "week", "month"), 
                                   units = -1, showRequestURL = FALSE) {
  unit_matched <- match.arg(unit)
  
  link_metrics_referrers_url <- "https://api-ssl.bitly.com/v3/link/referrers"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link,  
                limit = limit, unit = unit_matched, units = units)
  
  # call method from ApiKey.R
  df_link_metrics_referrers <- doRequest("GET", link_metrics_referrers_url, service = "bitly", query, showURL = showRequestURL)
  
  df_link_metrics_referrers_data <- df_link_metrics_referrers$data$referrers
  
  # sapply(df_link_metrics_referrers_data, class)
  
  return(df_link_metrics_referrers_data)
}


#' @title Returns metrics about the pages referring click traffic to a single Bitlink, grouped 
#' by referring domain.
#'  
#' @description \url{http://dev.bitly.com/link_metrics.html#v3_link_referrers_by_domain}
#' 
#' @inheritParams bitly_LinksMetricsClicks
#' 
#' @return clicks - the number of clicks referred from this domain.
#' @return referrer - the URL referring clicks.
#' 
#' @examples
#' \dontrun{
#' bitly_token <- bitly_auth(key = "", secret = "")
#' bitly_LinksMetricsReferrersByDomain("http://bit.ly/DPetrov",unit="day",units=-1,limit = 100)
#' }
#' 
#' @export
bitly_LinksMetricsReferrersByDomain <- function(link, limit = 1000, unit = c("minute", "hour", "day", "week", "month"), 
                                             units = -1, showRequestURL = FALSE) {
  unit_matched <- match.arg(unit)
  
  link_metrics_referrers_by_domain_url <- "https://api-ssl.bitly.com/v3/link/referrers_by_domain"
  
  query <- list(access_token = bitly_token$credentials$access_token, link = link, limit = limit, unit = unit_matched, units = units)
  
  # call method from ApiKey.R
  df_link_metrics_referrers_by_domain <- doRequest("GET", link_metrics_referrers_by_domain_url, 
                                                   service = "bitly", query, showURL = showRequestURL)
  
  df_link_metrics_referrers_by_domain_data <- df_link_metrics_referrers_by_domain$data$referrers
  
  # https://stackoverflow_com/questions/4227223/r-list-to-data-frame
  df_link_metrics_referrers_by_domain_data <- data.frame(t(sapply(df_link_metrics_referrers_by_domain_data, c)))
  
  df_link_metrics_referrers_by_domain_data$type <- rownames(df_link_metrics_referrers_by_domain_data) 
  df_link_metrics_referrers_by_domain_data$referrer <- as.character(df_link_metrics_referrers_by_domain_data$referrer)
  df_link_metrics_referrers_by_domain_data$clicks <- as.integer(df_link_metrics_referrers_by_domain_data$clicks)
  
  # sapply(df_link_metrics_referrers_by_domain_data, class)
  
  return(df_link_metrics_referrers_by_domain_data)
}
