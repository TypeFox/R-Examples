#' Authentication configuration/setup
#'
#' @export
#'
#' @family dsl
#' @name auth
#' @param .data Result of a call to \code{api}
#' @param token An OAuth token
#' @param app_name An OAuth application name
#' @param key An OAuth key
#' @param secret An OAuth secret key
#' @param base_url option url to use as base for request, authorize and access urls.
#' @param request url used to request initial (unauthenticated) token. If using
#' OAuth2.0, leave as NULL.
#' @param authorize url to send client to for authorisation
#' @param access url used to exchange unauthenticated for authenticated token.
#' @param user user name
#' @param pwd password
#' @param type type of HTTP authentication. Should be one of the following types
#' supported by Curl: basic, digest, digest_ie, gssnegotiate, ntlm, ntlm_vn, any.
#' Default: "basic" (the most common type)
#' @examples \dontrun{
#' # simple authentication (user/password)
#' api('https://httpbin.org/basic-auth/user/passwd') %>%
#'  api_simple_auth(user = "user", pwd = "passwd")
#' ## different auth type
#' # api('https://httpbin.org/basic-auth/user/passwd') %>%
#' #  api_simple_auth(user = "user", pwd = "passwd", type = "gssnegotiate")
#'
#' # OAuth setup
#' ## using a token
#' ### fill in your own token
#' # api('https://api.github.com/') %>%
#' #   api_oauth2(token = "<token>")
#'
#' # OAuth2
#' ## using a app name, key, and secret combination
#' ### uses a OAuth app set up by Hadley Wickham, which you'll auth against
#' # api('https://api.github.com/') %>%
#' #   api_oauth2(app_name = "github", key = "56b637a5baffac62cad9",
#' #       secret = "8e107541ae1791259e9987d544ca568633da2ebf",
#' #       base_url = "https://github.com/login/oauth",
#' #       authorize = "authorize", access = "access_token")
#'
#' # OAuth1
#' # api('https://api.twitter.com/1.1/statuses/home_timeline.json') %>%
#' #  api_oauth1(app_name = "twitter", key = "TYrWFPkFAkn4G5BbkWINYw",
#' #      secret = "qjOkmKYU9kWfUFWmekJuu5tztE9aEfLbt26WlhZL8",
#' #      base_url = "https://api.twitter.com/oauth/",
#' #      request = "request_token", authorize = "authenticate", access = "access_token")
#'
#' # Request some data with oauth2 via Github
#' ## put in your username and password
#' # api('https://api.github.com/') %>%
#' #  api_simple_auth(user = "<foo>", pwd = "<bar>")
#' }

# simple authentication ------------------------------------
#' @export
#' @rdname auth
api_simple_auth <- function(.data, user, pwd, type = "basic") {
  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  modifyList(.data, list(config = c(authenticate(user = user, password = pwd, type = type))))
}

# oauth ------------------------------------
# oauth2 ------------------------------------
#' @export
#' @rdname auth
api_oauth2 <- function(.data, token = NULL, app_name = NULL, key = NULL,
                       secret = NULL, base_url = NULL,
                       authorize = NULL, access = NULL) {

  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  args <- comp(list(token = token, app_name = app_name, key = key, secret = secret))
  if (length(args) == 0) {
    stop("either token or app_name + key + secret must be provided", call. = FALSE)
  } else {
    if (!is.null(token)) {
      auth <- config(token = token)
    } else {
      app <- oauth_app(app_name, key, secret)
      endpts <- oauth_endpoint(authorize = authorize, access = access, base_url = base_url)
      token <- oauth2.0_token(endpts, app)
      auth <- config(token = token)
    }
  }
  modifyList(.data, list(config = c(auth)))
}

# oauth1 ------------------------------------
#' @export
#' @rdname auth
api_oauth1 <- function(.data, token = NULL, app_name = NULL, key = NULL,
                       secret = NULL, base_url = NULL, request = NULL,
                       authorize = NULL, access = NULL) {

  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  args <- comp(list(token = token, app_name = app_name, key = key, secret = secret))
  if (length(args) == 0) {
    stop("either token or app_name + key + secret must be provided", call. = FALSE)
  } else {
    if (!is.null(token)) {
      auth <- config(token = token)
    } else {
      app <- oauth_app(app_name, key, secret)
      endpts <- oauth_endpoint(request = request, authorize = authorize,
                               access = access, base_url = base_url)
      token <- oauth1.0_token(endpts, app)
      auth <- config(token = token)
    }
  }
  modifyList(.data, list(config = c(auth)))
}
