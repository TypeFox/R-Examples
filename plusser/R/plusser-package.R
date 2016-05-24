#' plusser: A Google+ API Interface for R
#' 
#' The idea of this package is to provide a high level interface to Google's
#' Google+ social network. As of now, functions related to data retrieval are
#' available. The authentication (beyond using an API key) and posting
#' messages parts of the API are not yet implemented.
#' 
#' As a high-level interface between R and Google+, this package is aimed at
#' social media and researchers. Google does not provide an API to retrieve the
#' following relation ships without full authentication, this package is not yet
#' suitable for social network analysis.
#'
#' On the social media side, this package can be used to retrieve posts, a
#' post's popularity (comments, +1s, reshares) and the profiles of entities.
#' Espescially interesting might be the retrieval of profiles that interacted
#' with certain posts.
#' 
#' This package is designed as modular as possible with separating harvest and
#' parsing functions. Users can use their own parsing functions if required.
#'
#' @section Google+ API Key:
#'
#' This section describes briefly how to obtain a Google+ API key.
#' \enumerate{
#'   \item Go to the Google Developers Console at
#'     \url{https://console.developers.google.com/project}.
#'   \item Create a new project and open it.
#'   \item In the menu, choose \code{APIs & auth}.
#'   \item Choose Google+ API on the next screen.
#'   \item Activate it by clicking \code{On}.
#'   \item Choose credentials from the submenu.
#'   \item Click \code{Create new key} and write it down. This is your API key.
#' }
#'
#' @section Google+ Quotas:
#'
#' Currently, Google permits 5 requests per second up to a maximum total of
#' 10,000 requests per day.
#' @references See the official Google Google+ API documentation:
#'   \url{https://developers.google.com/+/api/latest/}.
#'
#' @session Suggested Workflow:
#'
#' The basic workflow is as follows:
#' \enumerate{
#'   \item Set you API key using \code{\link{setAPIkey}}.
#'   \item Locate the Google+ user ID or +name of the page you want to look up.
#'     Typically, by searching for that entity on Google+ and noting the UID or
#'     name.
#'   \item Harvest the page of that entitiy using \code{\link{harvestPage}}.
#'   \item Retrieve the UIDs of people resharing or +1ing a post with
#'     \code{\link{harvestActivity}}.
#'   \item Look up profile information for people acting on a post using
#'     \code{\link{harvestProfile}}.
#'   \item Use the resulting data frame as you see fit.
#' }
#' 
#' @import RCurl RJSONIO lubridate
#' @name plusser
#' @docType package
NULL
