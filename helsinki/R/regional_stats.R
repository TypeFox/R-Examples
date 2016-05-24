# This file is a part of the helsinki package (http://github.com/rOpenGov/helsinki)
# in association with the rOpenGov project (ropengov.github.io)

# Copyright (C) 2010-2014 Juuso Parkkinen / Louhos <louhos.github.com>. 
# All rights reserved.

# This program is open source software; you can redistribute it and/or modify 
# it under the terms of the FreeBSD License (keep this notice): 
# http://en.wikipedia.org/wiki/BSD_licenses

# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


#' Helsinki Region Infoshare statistics API
#'
#' Retrieves data from the Helsinki Region Infoshare (HRI) statistics API:
#' http://dev.hel.fi/stats/.
#' Currently provides access to the 'aluesarja't data: http://www.aluesarjat.fi/.
#' 
#' Current implementation is very simple.
#' You can either get the list of resources with query="",
#' or query for a specific resources and retrieve it in a
#' three-dimensional array form.
#'
#' @param query A string, specifying the dataset to query
#' @param verbose logical. Should R report extra information on progress? 
#' @return results A three-dimensional array of results.
#' 
#' @importFrom RCurl getCurlHandle
#' @importFrom RCurl getForm
#' @importFrom rjson fromJSON
#' @export
#' 
#' @references See citation("helsinki") 
#' @author Juuso Parkkinen \email{louhos@@googlegroups.com}
#' @examples stats.array <- get_hri_stats("aluesarjat_a03s_hki_vakiluku_aidinkieli")

get_hri_stats <- function (query="", verbose=TRUE) {
  
  ## TODO
  # implement grepping for resources? as in eurostat
  if (verbose)
    message("Accessing Helsinki Region Infoshare statistics API...")
  
  # Use the regional statistics API
  api.url <- "http://dev.hel.fi/stats/resources/"
  # For resources list
  if (query=="")
    query.url <- paste0(api.url, query)
  # For a specific resource, use jsontstat
  else
    query.url <- paste0(api.url, query, "/jsonstat")
  
  # Check whether url available
  if (!RCurl::url.exists(query.url)) {
    message(paste("Sorry! Url", query.url, "not available!\nReturned NULL."))
    return(NULL)
  }
  # Access data with RCurl
  curl <- RCurl::getCurlHandle(cookiefile = "")
  suppressWarnings(
    res.json <- RCurl::getForm(uri=query.url, curl=curl)
  )
  # Process json into a list
  res.list <- rjson::fromJSON(res.json)
  
  # Process and show list of resources
  if (query=="") {
    resources <- names(res.list[["_embedded"]])
    names(resources) <- sapply(res.list[["_embedded"]], function(x) x$metadata$label)
    if (verbose)
      message("Retrieved list of available resources.")
    return(resources)
    
  } else {
    
    ## Process jsonstat results into an array
    # For info about jsontstat, see http://json-stat.org/format/
    # Possible R package ot use: https://github.com/ajschumacher/rjstat
    
    # Process dimensions metadata
    dims <- res.list$dataset$dimension$size
    names(dims) <- res.list$dataset$dimension$id
    dimnames <- lapply(res.list$dataset$dimension[3:(length(dims)+2)], function(x) {res=unlist(x$category$label); names(res)=NULL; res})
    
    # Construct an array
    
    # For special characters:
    # Merkintojen selitykset:
    #   .. (kaksi pistetta), tietoa ei ole saatu, se on liian epavarma ilmoitettavaksi tai se on salattu;
    # . (piste), loogisesti mahdoton esitettavaksi;
    # 0 (nolla), suure pienempi kuin puolet kaytetysta yksikosta.
    # assign NA to ".", and ".."
    # => simple as.numeric() is fine, produces NA for "." and ".."
    
    # Have to reverse the dimensions, because in arrays
    # "The values in data are taken to be those in the array with the leftmost subscript moving fastest."
    res.list$dataset$value[res.list$dataset$value %in% c(".", "..")] <- NA
    res.array <- array(data=as.numeric(res.list$dataset$value), dim=rev(dims), dimnames=rev(dimnames))
    if (verbose)
      message("Retrieved resource '",query,"'")
    return(res.array)
  }
  #   # Test that it works
  #   query <- "aluesarjat_a03s_hki_vakiluku_aidinkieli"
  #   hki.vakiluku <- get_hri_stats(query)
  #   library(reshape2)
  #   df <- reshape2::melt(hki.vakiluku)
  
}
