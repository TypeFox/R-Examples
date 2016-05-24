#############################################################################
# 
# Copyright 2015 Axibase Corporation or its affiliates. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License.
# A copy of the License is located at
#
# https://www.axibase.com/atsd/axibase-apache-2.0.pdf
#
# or in the "license" file accompanying this file. This file is distributed
# on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
# express or implied. See the License for the specific language governing
# permissions and limitations under the License.
#
#############################################################################
#
#' Get information about entities from Axibase Time-Series Database.
#'
#' @description
#' This function fetches a list of entities from ATSD,
#' and convert it to a data frame.
#' 
#' @param expression
#'     Optional string argument.
#'     Select entities matching particular name pattern and/or user-defined entity tags.
#'     The syntax of the \code{expression} argument is explained in the package vignette.
#'     Type \code{browseVignettes(package = "atsd")} to see the vignette.
#' @param active
#'     Optional string argument: "true" or "false". 
#'     Filter entities by \code{lastInsertTime}. If 
#'     \code{active = "true"}, 
#'     only entities with positive lastInsertTime are included in the response.
#' @param tags
#'     Optional string argument.
#'     User-defined entitiy tags to be included in the response.
#'     By default, all the tags will be included.
#' @param limit
#'     Optional integer argument.
#'     If limit > 0, the response shows the top-N entities ordered by name.
#'     
#' @inheritParams query
#' 
#' @return 
#'     A data frame. 
#'     Each row of the data frame corresponds to an entity and its tags:
#'     \code{name}, 
#'     \code{enabled},
#'     \code{lastInsertTime}
#'     and user-defined entity tags as requested by the 'tags' argument.
#'     For more information look at the package vignette: 
#'     \code{browseVignettes(package = "atsd")}.
#' @seealso
#'     Visit
#'     \url{http://axibase.com/axibase-time-series-database/}
#'     for information about ATSD.
#' @examples \dontrun{
#' # get all entities and include all their tags in the data frame
#' get_entities()
#' 
#' # get all active entities and include all their tags in the data frame
#' get_entities(active = "true")
#' 
#' # Get the top 2 entities whose 'name'  and user-defined entity tag, 'app',
#' # match to the expression. Include the tag, 'app', into response 
#' # and exclude oter user-defined entity tags.
#' get_entities(expression = "name like 'nur*' and lower(tags.app) like '*hbase*'", 
#'              tags = "app", limit = 2)
#' }
#' @export
get_entities <- function(expression = "",
                        active = "",
                        tags = "*",
                        limit = "",
                        verbose = TRUE){
  
  if (check_connection()) {
    https_options <- set_https_options()
    response <- httr::GET(paste0(get("url", envir = atsdEnv), "/api/v1/entities"),
                          httr::authenticate(get("user", envir = atsdEnv), 
                                             get("password", envir = atsdEnv)),
                          query = list(
                            expression = expression,
                            active = active,
                            tags = tags,
                            limit = limit),
                          config = https_options
    )
    httr::stop_for_status(response)
    if (verbose) {
      cat("Your request was successfully processed by server. Start parsing and filtering.")
    }
    entities <- lapply(httr::content(response, "parsed"), make_flat)
    if (verbose) {
      cat("\nParsing and filtereng done. Start converting to data frame.")
    }
    
    # Conversion of list of entities to data frame.
    # The following one-liner is too slow, so we should do more work.
    # entities <- Reduce(merge_all, entities)
    entities <- list_to_dfr(entities)
    if (verbose) {
      cat("\nConverting to data frame done.")
    }
    
    if ("lastInsertTime" %in% names(entities)) {
      entities$lastInsertTime <- as.POSIXlt(entities$lastInsertTime / 1000, origin = "1970-01-01 00:00:00", tz = "GMT")
    }
    return(entities)
  } else {
    invisible()
  }
}
