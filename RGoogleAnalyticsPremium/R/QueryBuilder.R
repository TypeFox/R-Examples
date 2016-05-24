#' Initialize a QueryBuilder object with the given parameters and perform validation
#'
#' @export
#'
#' @param query.params.list List of all the Query Parameters. See \code{\link{Init}} for the
#' entire list
#'
#' @importFrom lubridate ymd
#'
#' @return The builder object to process the query parameters.

QueryBuilder <- function(query.params.list) {


  # Constants.
  kMaxDimensions <- 4
  kMaxMetrics <- 10


  # Query parameters.
  end.date     <- query.params.list$end.date
  metrics      <- query.params.list$metrics
  start.date   <- query.params.list$start.date
  title        <- query.params.list$title
  dimensions   <- query.params.list$dimensions
  filters      <- query.params.list$filters
  segment      <- query.params.list$segment


  #' Sets the Query Parameters for the Query Builder Object and performs validation
  #'
  #'

  SetQueryParams <- function() {

    EndDate(end.date)
    Metrics(metrics)
    StartDate(start.date)
    Title(title)
    Dimensions(dimensions)
    Filters(filters)
    Segment(segment)
    Title(title)

    #Perform Validation
    Validate()
    return(invisible())
  }


  #' Sets the start date.
  #' Optional.
  #' All Analytics feed requests must specify a beginning and ending date
  #' range. If you do not indicate start- and end-date values for the
  #' request, the server returns a request error.
  #' Date values are in the form YYYY-MM-DD.
  #' The earliest valid start-date is 2005-01-01. There is no upper limit
  #' restriction for a start-date. However, setting a start-date that is
  #' too  far in the future will most likely return empty results.
  #'
  #' @keywords internal
  #' @param
  #'    start.date.param Optional. A start date of the form "YYYY-MM-DD"
  #'                      as a string. If NULL is used, the start.date
  #'                      parameter will be unset. If no parameter is
  #'                      specified, the current start.date value is
  #'                      returned.
  #'
  #'  @return None The start.date value if start.date.param is not set.
  #'  Un-set the parameter if the value NULL is used.

  StartDate <- function(start.date.param = NA) {

    if (is.null(start.date.param)) {
      start.date <<- NULL
      return(invisible())
    }

    # Returns the current dimension value if no parameter is used.
    if (is.na(start.date.param)) {
      return(start.date)
    }

    # Error handling.
    # Check the form of the start.date.param.
    if (is.na(as.Date(start.date.param, "%Y-%m-%d"))) {
      stop("A start date must be specified of the form YYYY-MM-DD")
    }

    start.date <<- start.date.param
    return(invisible())
  }
  #' Sets the end date.
  #' Optional.
  #' All Analytics feed requests must specify a beginning and ending date
  #' range. If you do not indicate start- and end-date values for the
  #' request, the server returns a request error.
  #' Date values are in the form YYYY-MM-DD.
  #' The earliest valid start-date is 2005-01-01. There is no upper limit
  #' restriction for a start-date. However, setting a start-date that is
  #' too far in the future will most likely return empty results.
  #'
  #' @keywords internal
  #'  @param
  #'    end.date.param An end date of the form 'YYYY-MM-DD'
  #'                    as a string. If NULL is used, the end.date.param
  #'                    parameter will be unset. If no parameter is specified,
  #'                    the current end.date value is returned.
  #'
  #'  @return None The end.date value if end.date.param is not set.
  #'

  EndDate <- function(end.date.param = NA) {


    # Un-set the parameter if the value NULL is used.
    if (is.null(end.date.param)) {
      end.date <<- NULL
      return(invisible())
    }

    # Returns the current dimension value if no parameter is used.
    if (is.na(end.date.param)) {
      return(end.date)
    }

    # Error handling.
    # Check the form of the end.date.param.
    if (is.na(as.Date(end.date.param, "%Y-%m-%d"))) {
      stop("An end date must be specified of the form YYYY-MM-DD")
    }

    end.date <<- end.date.param
    return(invisible())
  }

  #' Sets the dimensions.
  #' @keywords internal
  #' Optional.
  #' The dimensions parameter defines the primary data keys for your
  #' Analytics report, such as ga:browser or ga:city. Use dimensions to
  #' segment your web property metrics. For example, while you can ask for
  #' the total number of pageviews to your site, it might be more
  #' interesting to ask for the number of pageviews segmented by browser.
  #' In this case, you'll see the number of pageviews from Firefox,
  #' Internet Explorer, Chrome, and so forth.
  #
  #' When using dimensions in a feed request, be aware of the following
  #' constraints:
  #'   You can supply a maximum of 4 dimensions for any query.
  #'   You can not send a query comprised only of dimensions:
  #'     You must combine any requested dimension with at least one metric.
  #'     Any given dimension can be used with other dimensions or metrics,
  #'       but only where Valid Combinations apply for that dimension.
  #'
  #' @seealso
  #' http://code.google.com/apis/analytics/docs/gdata/
  #'        gdataReferenceDimensionsMetrics.html#validCombinations
  #'
  #' NOTE: This method does not check for invalid dimensions or combinations.
  #'
  #'  @param
  #'    dimensions.param A vector of up to 4 dimensions, either as
  #'                      a single string or a vector or strings, E.g.
  #'                      "ga:source,ga:medium" or c("ga:source", "ga:medium")
  #'                      If NULL is used, the dimensions parameter will be
  #'                      unset. If no parameter is specified, the current
  #'                      dimension value is returned.
  #'
  #'  @return None The dimensions value if dimensions.param is not set.
  #'

  Dimensions <- function(dimensions.param = NA) {

    # Un-set the parameter if the value NULL is used.
    if (is.null(dimensions.param)) {
      dimensions <<- NULL
      return(invisible())
    }

    # Returns the current dimension value if no parameter is used.
    if (is.na(dimensions.param[1])) {
      return(dimensions)
    }

    # Error handling.
    # Validate the dimensions input is a vector.
    if (!is.vector(dimensions.param)) {
      stop(paste("dimensions must be a vector of string variables"))
    }

    # Error handling.
    # Validate the length of the vector is no greater than the max number
    # of allowed dimensions.
    if (length(dimensions.param) > kMaxDimensions) {
      stop(paste("Google Analytics can only handle up to", kMaxDimensions,
                 "dimensions parameters"))
    }

    # Error handling.
    # Validate the vector is a character type.
    # This will not stop a vector like, c("2", "this")
    if (!is.character(dimensions.param)) {
      stop(paste("dimensions must be character, please refer to the",
                 "Google Analytics API documentation for more information"))
    }

    dimensions <<- paste(dimensions.param, collapse = ",")
    return(invisible())
  }

  #' Sets the metrics of interest (clicks, pageviews, etc)
  #' Optional.
  #' The aggregated statistics for user activity in a profile, such as
  #' clicks or pageviews. When queried by alone, metrics provide aggregate
  #' values for the requested date range, such as overall pageviews or
  #' total bounces. However, when requested with dimensions, values are
  #' segmented by the dimension. For example, ga:pageviews requested with
  #' ga:country returns the total pageviews per country rather than the
  #' total pageviews for the entire profile. When requesting metrics, keep
  #' in mind:
  #'
  #' Any request must supply at least one metric because a request cannot
  #' consist only of dimensions.
  #' You can supply a maximum of 10 metrics for any query.
  #' Most combinations of metrics from multiple categories can be used
  #' together, provided no dimensions are specified.
  #' The exception to the above is the ga:visitors metric, which can only
  #' be used in combination with a subset of metrics.
  #' Any given metric can be used in combination with other dimensions or
  #' metrics, but only where Valid Combinations apply for that metric.
  #' Metric values are always reported as an aggregate because the Data
  #' Export API does not provide calculated metrics. For a list of common
  #' calculations based on aggregate metrics.
  #'
  #' NOTE: We do check for valid metrics.
  #'
  #' @keywords internal
  #' @param
  #'   metrics.param A vector of up to 10 metrics, either as
  #'                  a single string or a vector or strings. E.g.
  #'                  "ga:sessions" or c("ga:sessions", "ga:bounces")
  #'                  If NULL is used, the metrics parameter will be
  #'                  unset. If no parameter is specified, the current
  #'                  metrics value is returned.
  #' @return
  #'   The metrics value if metrics.param is not set.

  Metrics <- function(metrics.param = NA) {

    # Un-set the parameter if the value NULL is used.
    if (is.null(metrics.param)) {
      metrics <<- NULL
      return(invisible())
    }

    # Returns the current metrics value if no parameter is used.
    if (is.na(metrics.param[1])) {
      return(metrics)
    }

    # Error handling.
    # Check the metrics input is that of a vector
    if (!is.vector(metrics.param)) {
      stop("metrics must be a vector of string variables")
    }

    # Error handling.
    # Check the length of the vector is no greater than the max number of
    # metrics.
    if (length(metrics.param) > kMaxMetrics) {
      stop(paste("Google Analytics can only handle up to", kMaxMetrics,
                 "metrics parameters"))
    }

    # Error handling.
    # Check the vector is a character type.
    # this will not stop a vector like, c("2", "this")
    if (!is.character(metrics.param)) {
      stop(paste("metrics must be character string, please refer to the",
                 "Google Analytics API documentation for more information"))
    }

    # Combine and store the parameters.
    metrics <<- paste(metrics.param, collapse = ",")
    return(invisible())
  }

  #' Sets the segments, see dxp:segment in the Account Feed Response section
  #' in the GA literature online.
  #' http://code.google.com/apis/analytics/docs/gdata/gdataDeveloperGuide.html
  #' @keywords internal
  #' Optional.
  #' For general information on advanced segments, see Advanced
  #' Segmentation in the Help Center. You can request an advanced segment
  #' in the data feed in two ways:
  #'
  #' (1) The numeric ID of a default or custom advanced segment.
  #'     The account feed returns all default advanced segments and their
  #'     IDs, as well as any custom segments defined for the account.
  #'     For more information on segment and their IDs, see dxp:segment in
  #'     the Account Feed Response section.
  #' (2) The dynamic parameter in the query.
  #'     Use this method to segment your data request by one or more
  #'     dimensions and/or metrics. You can also use regular expressions
  #'     for segments just as you would for the filters parameter.
  #'     Dynamic segments use the same Expressions and Operators used for
  #'     the filters parameter. When using OR boolean logic or AND boolean
  #'     logic, dynamic segment expressions follow the same rules as for
  #'     the filters parameter, except that you may use OR boolean logic
  #'     with both dimensions or metrics.
  #' Dimensions/metrics combinations in the advanced segment expression
  #' have fewer restrictions. Except where noted in the table, you can use
  #' any dimension or metric in combination with another in your filter.
  #'
  #' The segment parameter is once again difficult to write checks for,
  #' as this is a handler we rely on the GA API to report errors with the
  #' request.
  #'
  #' Example:
  #' gaid::10
  #' dynamic::ga:medium==referral
  #' @param
  #'   segment An advanced segment definition to slice and dice your
  #'            Analytics data. If NULL is used, the segment parameter will be
  #'            unset. If no parameter is specified, the current segment value
  #'            is returned.
  #'
  #' @return None The segment value if segment.param is not set.
  #'

  Segment <- function(segment.param = NA) {

    # Un-set the parameter if the value NULL is used.
    if (is.null(segment.param)) {
      segment <<- NULL
      return(invisible())
    }

    # Returns the current segment value if no parameter is used.
    if (is.na(segment.param[1])) {
      return(segment)
    }

    segment <<- segment.param
    return(invisible())
  }

#   #' Sets the filters used.
  #' Optional.
  #' The filters query string parameter restricts the data returned from
  #' your request to the Analytics servers. When you use the filters
  #' parameter, you supply a dimension or metric you want to filter,
  #' followed by the filter expression. For example, the following feed
  #' query requests ga:pageviews and ga:browser from profile 12134, where
  #' the ga:browser dimension starts with the string Firefox:
  #'
  #' @param
  #'   filters The filter string for the GA request.
  #'            e.g. "ga:medium==referral".
  #'            If NULL is used, the filters parameter will be unset.
  #'            If no parameter is specified, the current filters value
  #'            is returned.
  #'
  #' @return None The filters value if filters.param is not set.
  #

  Filters <- function(filters.param = NA) {

    # Un-set the parameter if the value NULL is used.
    if (is.null(filters.param)) {
      filters <<- NULL
      return(invisible())
    }

    # Returns the current sort value if no parameter is used.
    if (is.na(filters.param[1])) {
      return(filters)
    }

    filters <<- filters.param
    return(invisible())
  }

#
  #' Returns whether the Query has all the required parameters set. These are
  #' the start.date, end.date, metrics, and table.id parameters.
  #'
  #' @return
  #'   TRUE if the query has all the required parameters. Otherwise stops the
  #'   program execution.
  Validate <- function() {

    missing.params <- c()

    if (is.null(start.date)) {
      missing.params <- append(missing.params, "start.date")
    }
    if (is.null(end.date)) {
      missing.params <- append(missing.params, "end.date")
    }
    if (is.null(metrics)) {
      missing.params <- append(missing.params, "metrics")
    }
    if (is.null(title)) {
      missing.params <- append(missing.params, "title")
    }
    if (length(missing.params) == 0) {
      return(TRUE)
    }
    missing.string <- paste(missing.params, collapse = ", ")
    stop(paste("All GA queries must have", missing.string, "parameters.",
               sep = " "))

    #Check whether start date precedes end date
    if ((as.numeric(difftime(ymd(start.date),ymd(end.date),units="days"))) < 0) {
      stop("Start Date must precede End date. Please correct the dates and re-initialize the query")
    }

  }



  #' A function to reset all the data values to NULL, for a new query.f
  #' The ClearData() function allows a user to reset the query parameters,
  #' (start.date, metrics, etc) back to NULL.
  #'@keywords internal
  #'
  #' @return
  #'  Resets all the query parameters to NULL.
  ClearData <- function() {

    start.date  <<- NULL
    end.date    <<- NULL
    dimensions  <<- NULL
    metrics     <<- NULL
    segment     <<- NULL
    filters     <<- NULL
    title       <<- NULL
    return(invisible())
  }

  #' Checks whether a valid title exists or not.
  #' @keywords internal
  #' @return None A stop call if the title is not valid or not present.
  #'
  Title <- function(title.param = NA) {

    if (is.null(title.param)) {
      title <<- NULL
      return(invisible())
    }

    # Returns the current value if no parameter is used.
    if (is.na(title.param[1])) {
      return(title)
    }

    # Error Handling.
    # Title must be string.
    if (!is.character(title.param)) {
      stop("title must be a string")
    }


    title <<- title.param

    return(invisible())
  }


  #' This function updates the start date in the query builder object
  #' @keywords internal
  #' @param Start Date
  #'
  #' @return None Sets the Start Date
  SetStartDate <- function(start.date) {
    StartDate(start.date)
  }

  GetStartDate <- function() {
    return(start.date)
  }

  #' This function updates the end date in the query builder object
  #' @param End Date
  #'
  #' @return None Sets the End Date
  SetEndDate <- function(end.date) {
    EndDate(end.date)
  }

  GetEndDate <- function() {
    return(end.date)
  }

  #' This function updates the start index in the query builder object
  #' @param End Date
  #'
  #' @return None Sets the Start Index in the Query Builder Object
#   SetStartIndex <- function(start.index) {
#     StartIndex(start.index)
#   }



  return(list("end.date"     =   EndDate,
              "metrics"      =   Metrics,
              "start.date"   =   StartDate,
              "title"        =   Title,
              "dimensions"   =   Dimensions,
              "filters"      =   Filters,
              "segments"     =   Segment,
              "clear.data"   =   ClearData,
              "Validate"     =   Validate,
              "GetStartDate"  = GetStartDate,
              "GetEndDate" = GetEndDate,
              "SetEndDate" = SetEndDate,
              "SetStartDate" = SetStartDate,
              "SetQueryParams" = SetQueryParams))
}

