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
#' Fetch time-series historic data or forecasts from Axibase Time-Series Database.
#'
#' @description
#' This function fetches time-series from ATSD and creates a data frame from it.
#' @param metric
#'     Required string argument.
#'     The name of the metric you want to get data for. 
#'     For example, \code{metric = "disk_used_percent"}.
#'     \cr
#'     To obtain a list of metrics collected by ATSD use the
#'     \code{\link{get_metrics}} function.
#' @param entity
#'     Optional string argument. 
#'     The name of the entity you want to get data for.
#'     If not provided, then data for all entities will be fetched
#'     for the specified metric.
#'     Obtain the list of entities and their tags with the
#'     \code{\link{get_entities}} function.
#' @param entity_group
#'     Optional string argument. 
#'     You could specify a group of entities and extract data 
#'     for entities from this group.
#'     For example, \code{entity_group = "HP Servers"}.
#' @param tags
#'     Optional string vector argument. 
#'     List of user-defined series tags to filter the fetched time-series data, 
#'     for example, 
#'     \code{c("disk_name=sda1", "mount_point=/")}.
#' @param selection_interval
#'     Required string argument.
#'     This is the time interval for which the data will be selected.
#'     Specify it as "n-unit", where  
#'     unit is a Second, Minute, Hour, Day, Week, Month, Quarter, or Year
#'     and n is the number of units, for example, "3-Week" or "12-Hour".
#' @param end_time
#'     Optional string argument.
#'     The end time of the selection interval, for example, 
#'     \code{end_time = "date(2014-12-27)"}.
#'     If not provided, the current time will be used. 
#'     Specify the date and time, or use one of the supported expressions:
#'     \url{http://axibase.com/axibase-time-series-database/screenshots-4/end-time/}.
#'     For example, \code{'current_day'} would set the end of selection interval 
#'     to 00:00:00 of the current day. 
#' @param aggregate_interval
#'     Optional string argument.
#'     The length of the aggregation interval. 
#'     The period of produced time-series will be equal to the 
#'     \code{aggregate_interval}.
#'     The value for each period is computed by the 
#'     \code{aggregate_statistics}
#'     function applied to all samples of the original time-series within the period.
#'     The format of the 
#'     \code{aggregate_interval} is the same as for the 
#'     \code{selection_interval} argument, for example, "1-Minute".
#' @param aggregate_statistics
#'     Optional string vector argument.
#'     The statistic function used for aggregation.
#'     List of available functions:
#'     "Avg", "Min", "Max", "Sum", "Count", "StDev", "WAvg", "WTAvg", 
#'     "Percentile 50", "Percentile 75", "Percentile 90", "Percentile 95", 
#'     "Percentile 99", "Percentile 99.5", "Percentile 99.9".
#'     Multiple values are supported, for example, c("Min", "Avg", "StDev").
#'     The default value is "Avg".
#' @param interpolation
#'     Optional string argument
#'     If aggregation is enabled, then the values for the periods without data 
#'     will be computed by one of the following interpolation functions: 
#'     "None", "Linear", "Step". The default value is "None".
#' @param export_type
#'     Optional string argument.
#'     It can take one of two values: "History" or "Forecast". The default value is "History".
#'     For example, \code{export_type = "Forecast"}.
#' @param verbose
#'     Optional boolean argument.
#'     If \code{verbose = FALSE} then all console output will be suppressed.
#'     By default, \code{verbose = TRUE}.
#' @return 
#'     The function returns a data frame. It could be empty if no data match your query
#'     or if your request could not be processed by ATSD server. In any case you will 
#'     get a console diagnostic message with a short description of the server responce.
#' @details
#'     The function has only two required arguments:
#'     \code{metric} and \code{selection_interval}.
#'     \cr
#'     Type 
#'     \code{browseVignettes(package = "atsd")}
#'     to view the complete package documentation and usage examples.
#' @seealso
#'     Visit
#'     \url{http://axibase.com/axibase-time-series-database/}
#'     for information about ATSD.
#' @examples \dontrun{
#' # Create data frame which contains time series for the given metric 
#' # and all entities for the last 1 hour.
#' dfr <- query(metric = "disk_used_percent", selection_interval = "1-Hour")
#' 
#' dfr <- query( export_type = "Forecast",
#'               metric = "disk_used_percent",
#'               entity_group = "Linux",
#'               tags = c("mount_point=/boot", "file_system=/dev/sda1"),
#'               selection_interval = "1-Week",
#'               aggregate_statistics = c("Avg", "Min", "Max"),
#'               aggregate_interval = "1-Minute",
#'               interpolation = "Linear")
#'                
#' # Example of the end_time argument usage.
#' dfr <- query( metric = "cpu_usage",
#'               entity = "host-383",
#'               selection_interval = "1-Day",
#'               end_time = "date('2015-02-10 10:15:03')")
#' }
#' @export
query <- function(metric,
                  entity = NA,
                  entity_group = NA,
                  tags = character(),
                  selection_interval,
                  end_time = NA,
                  aggregate_interval = NA,
                  aggregate_statistics = "Avg",
                  interpolation = "None",
                  export_type = "History",
                  verbose = TRUE) {

  missing_args <- check_arguments(export_type, metric, selection_interval)
  if (!check_connection() || missing_args != "") {
    if (verbose) {
      cat(missing_args)
    }
    return()
  }
  req <- get_request( export_type = export_type,
                      metric = metric,
                      entity = entity,
                      entity_group = entity_group,
                      tags = tags,
                      selection_interval = selection_interval,
                      end_time = end_time,
                      aggregate_interval = aggregate_interval,
                      interpolation = interpolation,
                      aggregate_statistics = aggregate_statistics)
  request <- req[1]
  warnings <- req[2]
  if (verbose) {
    cat(warnings)
  }
  
  url <- get("url", envir = atsdEnv)
  userpwd <- paste0(get("user", envir = atsdEnv), ":", get("password", envir = atsdEnv))
  url_encoded <- RCurl::curlEscape(request)
  request <- paste0(url, '/export?settings=', request)
  url_encoded <- paste0(url, '/export?settings=', url_encoded)
  
  
  d <- RCurl::debugGatherer()
  
  https_options <- set_https_options()
  
  response_body <- RCurl::getURL(url_encoded,
                                 userpwd = userpwd, 
                                 verbose = TRUE, 
                                 httpauth = 1L,
                                 debugfunction = d$update,
                                 .opts = https_options)
  
  response_header <- unname(d$value()["headerIn"])
  response <- parse_response(response_body, response_header)
  if (verbose) {
    cat(response[[2]])
  }
  return(response[[1]])
}
