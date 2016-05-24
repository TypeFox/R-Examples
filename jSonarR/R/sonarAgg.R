#
#  Copyright 2014 jSonar Inc
#  All Rights Reserved.
#
#  Licensed under the GNU AFFERO GENERAL PUBLIC LICENSE version 3
#  See http://www.r-project.org/Licenses/AGPL-3
#

#' Run a saved aggregation pipeline
#'
#' Execute an aggregation pipeline which has been saved and published in JSON
#' Studio Analytics, and get the result in a data frame.
#'
#' The parameters for this function are explained in greater detail in the
#' JSON Studio help page \emph{Using the Gateway}.
#'
#' @param connection a SonarConnection object created with
#'   \code{\link{new.SonarConnection}}
#' @param queryName the name of the saved query to execute
#' @param queryCol a collection in the database to use with the query
#' @param bind a list of bind variables and their values
#' @param limit the maximum number of results to return
#' @param idCol the name of a field which uniquely identifies each document.
#'   This will be used for the row names in the returned data frame. The
#'   default is X_id, which is the name of Mongo's _id field (adjusted by
#'   \code{\link{make.names}}).
#' @param publishedBy the name of the user who we expect published the API
#' @param colClasses a list of column names and their respective classes, as
#'   used in \code{\link{read.csv}}. This may be necessary if some columns'
#'   types are not being detected automatically.
#'
#' @examples
#' connection <- new.SonarConnection('https://example.com', 'localhost', 'test')
#'
#' ny_by_day <- sonarAgg(connection, 'delays_by_day', 'NYCFlights')
#' cor(ny_by_day$X_avg_ArrDelay, ny_by_day$X_avg_AirTime) 
#'
#' @export
#' @keywords database
#'
#' @family connection
#' @family csv
#' @seealso \url{http://jsonstudio.com/wp-content/uploads/2014/04/manual141/_build/html/index.html}
sonarAgg <- function(connection, queryName, queryCol, bind=list(), limit=NULL, idCol="_id", publishedBy=NULL, colClasses=NA)
{
    return(sonarCSV(connection, queryName, queryCol, type='agg', bind, limit, idCol, publishedBy, colClasses));
}
