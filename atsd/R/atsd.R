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
#' Support querying Axibase Time-Series Database.
#'
#' @description
#' The package lets you query 
#' \href{http://axibase.com/axibase-time-series-database/}{Axibase Time-Series Database}
#' (ATSD)
#' for time-series data and forecasts.
#' List of package functions:
#'
#' \code{\link{set_connection}, \link{save_connection}, 
#' \link{show_connection}} - are used to manage connection with ATSD: 
#' set up and store the url, user name, and password,
#' configure cryptographic protocol and enforce SSL certificate validation 
#' in the case of https connection.
#'
#' \code{\link{query}} -- get time-series data and forecasts from ATSD.
#'
#' \code{\link{get_metrics}} -- get information about the metrics collected by ATSD.
#' 
#' \code{\link{get_entities}} -- get information about the entities collected by ATSD.
#' 
#' \code{\link{to_zoo}} - converts a time-series data frame to 'zoo' object
#'  for manipulating irregular time-series with built-in functions in zoo package.
#' 
#' Type 
#' \strong{\code{browseVignettes(package = "atsd")}}
#' to view the complete package documentation and usage examples.
#' 
#' @author Axibase, api@@axibase.com
#' @docType package
#' @name atsd
NULL
