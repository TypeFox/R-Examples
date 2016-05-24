# The inner package environment
.RGAEnv <- new.env(parent = emptyenv())
.RGAEnv$Token <- NULL

# Set token to environment
set_token <- function(value) {
    .RGAEnv$Token <- value
    return(value)
}

# Get token from environment
get_token <- function() {
    .RGAEnv$Token
}

# Remove token from environment
remove_token <- function(token) {
    cache_path <- token$cache_path
    default_path <- getOption("rga.cache")
    if (is.null(cache_path))
        cache_path <- default_path
    if (!is.null(cache_path) && cache_path != default_path)
        cache_path <- c(cache_path, default_path)
    e <- file.exists(cache_path)
    if (any(e)) {
        message(sprintf("Removed old cache files: %s.", paste(cache_path[e], collapse = ", ")))
        file.remove(cache_path[e])
    }
    set_token(NULL)
}

# Validate token
validate_token <- function(token) {
    if (missing(token))
        stop("Authorization error. Access token not found.")
    if (!inherits(token, "Token2.0"))
        stop(sprintf("Token is not a Token2.0 object. Found: %s.", class(token)))
    if(!token$validate())
        token$refresh()
    token$validate()
}

# Check environment variables exists
env_exists <- function(...) {
    dots <- list(...)
    nzchar(Sys.getenv(dots))
}

# Fix username without domain
fix_username <- function(x) {
    if (!grepl("@", x, fixed = TRUE))
        x <- paste0(x, "@gmail.com")
    return(x)
}

#' @title Authorize the RGA package to the user's Google Analytics account using OAuth2.0
#'
#' @description \code{authorize()} function uses \code{\link[httr]{oauth2.0_token}} to obtain the OAuth tokens. Expired tokens will be refreshed automamaticly. If you have no \code{client.id} and \code{client.secret} the package provides predefined values.
#'
#' @param username character. Google username email address hint. If not set you will need choose an account for the authorization.
#' @param client.id character. OAuth client ID. If you set the environment variable \env{RGA_CLIENT_ID} it is used.
#' @param client.secret character. OAuth client secret. If you set the environment variable \env{RGA_CLIENT_SECRET} it is used.
#' @param cache logical or character. \code{TRUE} means to cache using the default cache file \code{.oauth-httr}, \code{FALSE} means not to cache. A string means to use the specified path as the cache file. Otherwise will be used the \code{rga.cache} option value (\code{.ga-token.rds} by default). If \code{username} argument specified token will be cached in the \code{.username-token.rds} file.
#' @param reauth logical. Set \code{TRUE} to reauthorization with the same or different Google Analytics account.
#' @param token A valid \code{Token2.0} object (icluding \code{TokenServiceAccount}) to setup directly.
#'
#' @details
#'
#' After calling this function first time, a web browser will be opened. First, log in with a Google Account, confirm the authorization to access the Google Analytics data. Note that the package requests access for read-only data.
#'
#' When the \code{authorize()} function is used the \code{Token} variable is created in the separate \code{.RGAEnv} environment which is not visible for user. So, there is no need to pass the token argument to any function which requires authorization every time. Also there is a possibility to store token in separate variable and to pass it to the functions. It can be useful when you are working with several accounts at the same time.
#'
#' \code{username}, \code{client.id} and \code{client.secret} params can be specified by an appropriate options (with \dQuote{RGA} prefix): \env{RGA_USERNAME}, \env{RGA_CLIENT_ID}, \env{RGA_CLIENT_SECRET}.
#'
#' @section Use custom Client ID and Client secret:
#'
#' Google Analytics is used by millions of sites. To protect the system from receiving more data than it can handle, and to ensure an equitable distribution of system resources, certain limits have been put in place.
#'
#' The following quota limits are shared between all \pkg{RGA} users which use the predefined credentials (daily quotas refresh at midnight PST):
#'
#' \itemize{
#'   \item 50,000 requests per day
#'   \item 10 queries per second per IP
#' }
#'
#' To get full quota, you must register your application in the Google Developers Console. When you register a new application, you are given a unique client ID to identify each application under that project. To find your project's client ID and client secret, do the following:
#'
#' \enumerate{
#'   \item Open the \href{https://console.developers.google.com/projectselector/apis/credentials}{Credentials page}.
#'   \item Select a project (create if needed).
#'   \item create your project's OAuth 2.0 credentials by clicking \emph{Add credentials} > \emph{OAuth 2.0 client ID} and select \emph{Other} application type.
#'   \item Look for the Client ID in the OAuth 2.0 client IDs section. You can click the application name for details.
#' }
#'
#' To enable Analytics API for your project, do the following:
#'
#' \enumerate{
#'   \item Open the \href{https://console.developers.google.com/projectselector/apis/api/analytics/overview}{Analytics API Overview page}.
#'   \item CLick on the \emph{Enable API} button to activate Analytics API.
#' }
#'
#' There 3 ways to use custom Client ID and Client secret:
#'
#' \enumerate{
#'   \item Pass the \code{client.id} and \code{client.secret} arguments directly in the \code{authorize()} function call
#'   \item Set the \env{RGA_CLIENT_ID} and \env{RGA_CLIENT_SECRET} environment variables
#'   \item Set the \code{rga.client.id} and \code{rga.client.secret} options into the R session
#' }
#'
#' @section Revoke access application:
#'
#' To revoke access the \pkg{RGA} package do the following:
#'
#' \enumerate{
#'   \item Go to the \href{https://security.google.com/settings/security/permissions}{Apps connected to your account} page
#'   \item Find \emph{RGA package} entry. Then click on it
#'   \item Click on the \emph{Revoke access} button in the sidebar on the right
#' }
#'
#' @return A \code{\link[httr]{Token2.0}} object containing all the data required for OAuth access.
#'
#' @references \href{https://console.developers.google.com/}{Google Developers Console}
#'
#' \href{http://en.wikipedia.org/wiki/Environment_variable}{Environment variable}
#'
#' @seealso
#'
#' Other OAuth: \code{\link[httr]{oauth_app}} \code{\link[httr]{oauth2.0_token}} \code{\link[httr]{Token-class}}
#'
#' To revoke all tokens: \code{\link[httr]{revoke_all}}
#'
#' Setup environment variables: \code{\link{Startup}}
#'
#' @examples
#' \dontrun{
#' authorize(client.id = "my_id", client.secret = "my_secret")
#' # if set RGA_CLIENT_ID and RGA_CLIENT_SECRET environment variables
#' authorize()
#' # assign token to variable
#' ga_token <- authorize(client.id = "my_id", client.secret = "my_secret")
#' }
#'
#' @importFrom httr oauth_app oauth_endpoints oauth2.0_token
#' @export
authorize <- function(username = getOption("rga.username"),
                      client.id = getOption("rga.client.id"),
                      client.secret = getOption("rga.client.secret"),
                      cache = getOption("rga.cache"),
                      reauth = FALSE, token = NULL) {
    if (!is.null(token) && validate_token(token)) {
        set_token(token)
        return(invisible(token))
    }
    if (all(env_exists("RGA_CLIENT_ID", "RGA_CLIENT_SECRET"))) {
        message("Client ID and Client secret loaded from environment variables.")
        client.id <- Sys.getenv("RGA_CLIENT_ID")
        client.secret <- Sys.getenv("RGA_CLIENT_SECRET")
    }
    if (env_exists("RGA_USERNAME")) {
        message("Username loaded from environment variable.")
        username <- Sys.getenv("RGA_USERNAME")
    }
    app <- oauth_app(appname = "rga", key = client.id, secret = client.secret)
    endpoint <- oauth_endpoints("google")
    if (!is.null(username)) {
        stopifnot(is.character(username))
        stopifnot(length(username) == 1)
        username <- fix_username(username)
        endpoint$authorize <- paste0(endpoint$authorize, "?login_hint=", username)
        if (is.character(cache))
            cache <- paste0(".", username, "-token.rds")
    }
    if (reauth)
        remove_token(get_token())
    if (is.character(cache) && nzchar(cache))
        message(sprintf("Access token will be stored in the '%s' file.", cache))
    token <- oauth2.0_token(endpoint = endpoint, app = app, cache = cache,
                                  scope = "https://www.googleapis.com/auth/analytics.readonly")
    if (validate_token(token))
        set_token(token)
    return(invisible(token))
}
