get_auth_url <- function(creds) {
    url <- httr::modify_url(creds$auth_uri,
                            query = list(
                                client_id = creds$client_id,
                                scope = "https://www.googleapis.com/auth/analytics.readonly",
                                redirect_uri = creds$redirect_uris,
                                response_type = "code"))
    return(url)
}

get_auth_token <- function(creds, code) {
    endpoint <- httr::oauth_endpoints("google")
    req <- httr::POST(endpoint$access, encode = "form",
                      body = list(
                          client_id = creds$client_id,
                          client_secret = creds$client_secret,
                          redirect_uri = creds$redirect_uris,
                          grant_type = "authorization_code",
                          code = code))
    resp <- httr::content(req)
    app <- httr::oauth_app("ga-auth", creds$client_id, creds$client_secret)
    token <- httr::Token2.0$new(app = app, endpoint = endpoint,
                                credentials = resp, cache_path = FALSE,
                                params = list(scope = "https://www.googleapis.com/auth/analytics.readonly",
                                              as_header = TRUE))
    return(token)
}

creds <- jsonlite::fromJSON("creds.json")$web
auth_url <- get_auth_url(creds)
delay <- 10000 # in miliseconds

id <- 7783180
