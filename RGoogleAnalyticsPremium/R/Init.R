#' Initialize the Google Analytics query parameters
#' 
#' This function takes all the query parameters and combines them into a single
#' list that is to be passed as an argument to \code{\link{QueryBuilder}}. Note
#' that parameter validation is performed when the \code{\link{QueryBuilder}}
#' object is created
#' 
#' @export
#' @param start.date Start Date for fetching Analytics Data. Start Date must be
#'   of the format "\%Y-\%m-\%d"
#'   
#' @param end.date End Date for fetching Analytics Data. End Date must be of the
#'   format "\%Y-\%m-\%d"
#'   
#'   
#' @param metrics A vector of up to 10 metrics, either as a single string or a
#'   vector or strings. E.g. "ga:sessions" or c("ga:sessions", "ga:bounces").
#'
#' @param title Title of unsampled report.
#'
#' @param dimensions Optional. A vector of up to 4 dimensions, either as a
#'   single string or a vector or strings, E.g. "ga:source,ga:medium" or
#'   c("ga:source", "ga:medium").
#'   
#' @param filters Optional.The filter string for the GA request.e.g. "ga:medium==referral".
#' 
#' @param segments Optional.An advanced segment definition to slice and dice
#'   your Analytics data.
#'   
#' @return List of all the Query Parameters initialized by the user

Init <- function(
  end.date = NULL,
  metrics = NULL,
  start.date = NULL,
  title = NULL,
  dimensions = NULL,
  filters = NULL,
  segments = NULL){
  
  query.params.list = list("end.date" = end.date,
                           "metrics" = metrics,
                           "start.date" = start.date,
                           "title" = title,
                           "dimensions" = dimensions,
                           "filters" = filters,
                           "segments" = segments)
  
  return(query.params.list)
}
