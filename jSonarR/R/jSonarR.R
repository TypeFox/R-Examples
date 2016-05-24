#
#  Copyright 2014 jSonar Inc
#  All Rights Reserved.
#
#  Licensed under the GNU AFFERO GENERAL PUBLIC LICENSE version 3
#  See http://www.r-project.org/Licenses/AGPL-3
#

#' jSonar Analytics Platform API for R
#'
#' This package enables users to access MongoDB by running queries
#' and returning their results in R data frames. Usually, data in MongoDB is
#' only available in the form of a JSON document. jSonarR uses data
#' processing and conversion capabilities in the jSonar Analytics Platform
#' and the JSON Studio Gateway (\url{http://www.jsonstudio.com}), to convert it to
#' a tabular format which is easy to use with existing R packages.
#'
#' To use jSonarR, you must have access to a server running JSON Studio.
#' Create a connection using using \code{\link{new.SonarConnection}}. Now
#' you can run a saved query against a collection in the database using the
#' connection object and \code{\link{sonarAgg}} or \code{\link{sonarFind}}.
#'
#' @examples
#' connection <- new.SonarConnection('https://example.com', 'localhost', 'test')
#'
#' ny_by_day <- sonarAgg(connection, 'delays_by_day', 'NYCFlights')
#' summary(ny_by_day)
#'
#' tx_to_co <- sonarFind(connection, 'flights_to', 'TXFlights',
#'   bind=list(state="CO"),
#'   colClasses=c(DAY_OF_MONTH='factor', DEST_AIRPORT_ID='factor'))
#' summary(tx_to_co$DEST_AIRPORT_ID)
#'
#' @seealso MongoDB \url{http://www.mongodb.org}
#' @seealso JSON Studio \url{http://www.jsonstudio.com}
#' @import RCurl
#' @import methods
#' @import jsonlite
#' @docType package
#' @keywords connection database
#' @name jSonarR
NULL
