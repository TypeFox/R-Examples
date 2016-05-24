#' Discgolf setup
#'
#' @name disc_setup
#' @section How to setup:
#' You can set your url, username, and api key either using R options or
#' environment variables. I recommend using environment variables since
#' they are more general to any programming language, and you can easily
#' set secure env vars on a server or e.g., if you're doing continuous
#' integration on Travis-CI (or elsewhere).
#'
#' @section URL:
#' The base URL for the Discourse instance, e.g., \code{https://meta.discourse.org}
#'
#' Use:
#' \itemize{
#'  \item Env var name: \code{DISCOURSE_URL}
#'  \item R option name: \code{discourse_url}
#' }
#'
#' @section Username:
#' The user name you have registered on the Discourse instance you want
#' to use.
#'
#' Use:
#' \itemize{
#'  \item Env var name: \code{DISCOURSE_USERNAME}
#'  \item R option name: \code{discourse_username}
#' }
#'
#' @section API key:
#' The API key on the Discourse instance you want to use. This is not your
#' password you used to login to the instance. If you're the admin, you can get
#' an API key by going the dashboard at \code{base_url/admin}, then the API tab
#' \code{base_url/admin/api}, then generate a key, or copy the one already there.
#'
#' Use:
#' \itemize{
#'  \item Env var name: \code{DISCOURSE_API_KEY}
#'  \item R option name: \code{discourse_api_key}
#' }
NULL
