#' @title A Google Analytics API client for R
#'
#' @description
#' A package for extracting data from Google Analytics API into R.
#'
#' @section Key features:
#'
#' \itemize{
#'   \item Support for \href{https://developers.google.com/accounts/docs/OAuth2}{OAuth 2.0 authorization};
#'   \item Access to the following \href{https://developers.google.com/analytics/devguides/platform/}{Google Analytics APIs}:
#'   \itemize{
#'     \item \href{https://developers.google.com/analytics/devguides/config/mgmt/v3}{Management API}: access to configuration data for accounts, web properties, views (profiles), goals and segments;
#'     \item \href{https://developers.google.com/analytics/devguides/reporting/core/v3}{Core Reporting API}: query for dimensions and metrics to produce customized reports;
#'     \item \href{https://developers.google.com/analytics/devguides/reporting/mcf/v3}{Multi-Channel Funnels Reporting API}: query the traffic source paths that lead to a user's goal conversion;
#'     \item \href{https://developers.google.com/analytics/devguides/reporting/realtime/v3}{Real Time Reporting API}: report on activity occurring on your property at the moment;
#'     \item \href{https://developers.google.com/analytics/devguides/reporting/metadata/v3}{Metadata API}: access the list of API dimensions and metrics and their attributes;
#'   }
#'   \item Access to all the accounts which the user has access to;
#'   \item API responses is converted directly into R as a \code{data.frame};
#'   \item Auto-pagination to return more than 10,000 rows of the results by combining multiple data requests.
#' }
#'
#' To report a bug please type: \code{utils::bug.report(package = "RGA")}.
#'
#' @section Useage:
#'
#' Once you have the package loaded, there are 3 steps you need to use to get data from Google Analytics:
#'
#' \enumerate{
#'   \item Authorize this package to access your Google Analytics data with the \code{\link{authorize}} function;
#'   \item Determine the profile ID which you want to get access to with the \code{\link{list_profiles}} function;
#'   \item Get the results from the API with one of these functions: \code{\link{get_ga}}, \code{\link{get_mcf}} or \code{\link{get_realtime}}.
#' }
#'
#' For details about this steps please type into R: \code{browseVignettes(package = "RGA")}
#'
#' @section Bug reports:
#'
#' Before posting a bug please try execute your code with the \code{\link[httr]{with_verbose}} wrapper. It will be useful if you attach verbose output to the bug report. For example: \code{httr::with_verbose(list_profiles())}
#'
#' Post the \code{traceback()} output also may be helpful.
#'
#' To report a bug please type into R: \code{utils::bug.report(package = "RGA")}
#'
#' @name RGA-package
#' @docType package
#' @keywords package
#' @aliases rga RGA
#'
#' @author Artem Klevtsov \email{a.a.klevtsov@@gmail.com}
#'
#' @examples
#' \dontrun{
#' # load package
#' library(RGA)
#' # get access token
#' authorize()
#' # get a GA profiles
#' ga_profiles <- list_profiles()
#' # choose the profile ID by site URL
#' id <- ga_profiles[grep("http://example.com", ga_profiles$website.url), "id"]
#' # get date when GA tracking began
#' first.date <- firstdate(id)
#' # get GA report data
#' ga_data <- get_ga(id, start.date = first.date, end.date = "today",
#'                   metrics = "ga:users,ga:sessions",
#'                   dimensions = "ga:userGender,ga:userAgeBracket")
#' }
#'
NULL
