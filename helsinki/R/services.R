# This file is a part of the helsinki package (http://github.com/rOpenGov/helsinki)
# in association with the rOpenGov project (ropengov.github.io)

# Copyright (C) 2010-2014 Juuso Parkkinen, Leo Lahti and Joona Lehtomaki / Louhos <louhos.github.com>. 
# All rights reserved.

# This program is open source software; you can redistribute it and/or modify 
# it under the terms of the FreeBSD License (keep this notice): 
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Access Helsinki region Service Map API
#'
#' Access the new Helsinki region Service Map (Paakaupunkiseudun Palvelukartta)
#' http://dev.hel.fi/servicemap/ data through the API: http://api.hel.fi/servicemap/v1/. 
#' For more API documentation and license information see the API link.
#' 
#' @param query The API query as a string, for example "search", "service", or "unit".
#' For full list of available options and details, see http://api.hel.fi/servicemap/v1/. 
#' @param ... Additional parameters to the API (optional).
#' For details, see http://api.hel.fi/servicemap/v1/. 
#'
#' @return List of results
#' @export
#' @importFrom RCurl getCurlHandle
#' @importFrom RCurl getForm
#' @importFrom rjson fromJSON
#' 
#' @author Juuso Parkkinen \email{louhos@@googlegroups.com}
#' @examples search.puisto <- get_servicemap(query="search", q="puisto")

get_servicemap <- function(query, ...) {
  
  # api.url <- "http://www.hel.fi/palvelukarttaws/rest/v2/"
  # Define query url
  # New API (13.5.2014)
  api.url <- "http://api.hel.fi/servicemap/v1/"
  
  # Check whether API url available
  if (!RCurl::url.exists(api.url)) {
    message(paste("Sorry! API", api.url, "not available!\nReturned NULL."))
    return(NULL)
  }
  
  query.url <- paste0(api.url, query, "/")
  
  # Get Curl handle
  curl <- RCurl::getCurlHandle(cookiefile = "")
  
  # Get data as json using getForm
  # Note! Warnings suppressed because getForm outputs warning when no parameters (...) given
  suppressWarnings(
    res.json <- RCurl::getForm(uri=query.url, ..., curl=curl)
  )
  # Transform results into list from JSON
  res.list <- rjson::fromJSON(res.json)
  return(res.list)
}



#' Access Helsinki Linked Events API
#'
#' Access the new Helsinki Linked Events API: http://api.hel.fi/linkedevents/v0.1/.
#' The API contains data from the Helsinki City Tourist & Convention Bureau, 
#' the City of Helsinki Cultural Office and the Helmet metropolitan area public libraries.
#' For more API documentation and license information see the API link.
#' 
#' @param query The API query as a string, one of "category", "event", "language", or "place".
#' For details, see http://api.hel.fi/linkedevents/v0.1/. 
#' @param ... Additional parameters to the API (optional).
#' For details, see http://api.hel.fi/linkedevents/v0.1/. 
#'
#' @return List of results
#' @export
#' @importFrom RCurl getCurlHandle
#' @importFrom RCurl getForm
#' @importFrom rjson fromJSON
#' 
#' @author Juuso Parkkinen \email{louhos@@googlegroups.com}
#' @examples events <- get_linkedevents(query="event")

get_linkedevents <- function(query, ...) {
  
  # Define query url
  api.url <- "http://api.hel.fi/linkedevents/v0.1/"
  
  # Check whether API url available
  if (!RCurl::url.exists(api.url)) {
    message(paste("Sorry! API", api.url, "not available!\nReturned NULL."))
    return(NULL)
  }
  
  query.url <- paste0(api.url, query, "/")
  
  # Get Curl handle
  curl <- RCurl::getCurlHandle(cookiefile = "")
  
  # Get data as json using getForm
  # Note! Warnings suppressed because getForm outputs warning when no parameters (...) given
  suppressWarnings(
    res.json <- RCurl::getForm(uri=query.url, ..., curl=curl)
  )
  # Transform results into list from JSON
  res.list <- rjson::fromJSON(res.json)
  return(res.list)
}


