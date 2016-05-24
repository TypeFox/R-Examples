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
#' Get information about metrics from Axibase Time-Series Database.
#'
#' @description
#' This function fetches a list of metrics and their tags from ATSD,
#' and converts it to a data frame.
#' @param expression
#'     Optional string argument.
#'     Select metrics matching particular name pattern and/or user-defined metric tags.
#'     The syntax of the \code{expression} argument is explained in the package vignette.
#'     Type \code{browseVignettes(package = "atsd")} to see the vignette.
#' @param active
#'     Optional string argument: "true" or "false". 
#'     Filter metrics by \code{lastInsertTime}. If 
#'     \code{active = "true"}, 
#'     only metrics with positive lastInsertTime are included in the response.
#' @param tags
#'     Optional string argument.
#'     User-defined metric tags to be included in the response.
#'     By default, all the tags will be included.
#' @param limit
#'     Optional integer argument.
#'     If limit > 0, the response shows the top-N metrics ordered by name.
#' 
#' @inheritParams query
#' 
#' @return 
#'     A data frame. 
#'     Each row of the data frame corresponds to a metric and its tags:
#'     \code{name}, 
#'     \code{counter},
#'     \code{lastInsertTime}
#'     and user-defined metric tags as requested by the 'tags' argument.
#'     For more information view the package vignette: 
#'     \code{browseVignettes(package = "atsd")}.
#' @seealso
#'     Visit
#'     \url{http://axibase.com/axibase-time-series-database/}
#'     for information about ATSD.
#' @examples \dontrun{
#' # get all metrics and include all their tags in the data frame
#' get_metrics()
#' 
#' # get the top 100 active metrics which have tag, 'table', 
#' # include this tag into response and exclude oter user-defined metric tags
#' get_metrics(expression = "tags.table != ''", active = "true", 
#'             tags = "table", limit = 100)
#' 
#' # get metrics which have user-defined metric tag, 'table', 
#' # and whose name starts with 'cpu'
#' get_metrics(expression = "name like 'cpu*' and tags.table != ''")
#'             
#' # more complitcated expressions
#' get_metrics(expression = "likeAll(name, list('*disk*,*use*'))")
#' get_metrics(expression = "(name like 'cpu*' or tags.source = '') and tags.table like 'BC*'")
#' }
#' @export
get_metrics <- function(expression = "",
                        active = "",
                        tags = "*",
                        limit = "",
                        verbose = TRUE){
  
  if (check_connection()) {
    https_options <- set_https_options()
    response <- httr::GET(paste0(get("url", envir = atsdEnv), "/api/v1/metrics"),
                          httr::authenticate(get("user", envir = atsdEnv), 
                                             get("password", envir = atsdEnv)),
                          query = list(
                            expression = expression,
                            active = active,
                            tags = tags,
                            limit = limit),
                          config = https_options
                            #list(
                            #cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl"))), 
                            #followlocation = TRUE, 
                            #useragent = "R", 
                            #ssl.verifypeer = FALSE
                            #sslversion=3)
    )
    httr::stop_for_status(response)
    if (verbose) {
      cat("Your request was successfully processed by server. Start parsing and filtering.")
    }
    metrics <- lapply(httr::content(response, "parsed"), filter_metric_attributes)
    metrics <- lapply(metrics, make_flat)
    if (verbose) {
      cat("\nParsing and filtereng done. Start converting to data frame.")
    }
    
    # Conversion of list of metrics to data frame.
    # The following one-liner is too slow, so we should do more work.
    # metrics <- (Reduce(merge_all, metrics))
    metrics <- list_to_dfr(metrics)
    if (verbose) {
      cat("\nConverting to data frame done.")
    }
    
    if ("lastInsertTime" %in% names(metrics)) {
      metrics$lastInsertTime <- as.POSIXlt(metrics$lastInsertTime / 1000, origin = "1970-01-01 00:00:00", tz = "GMT")
    }
    return(metrics)
  } else {
    invisible()
  }
}
