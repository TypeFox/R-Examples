.state <- new.env(parent = emptyenv())
globalVariables(c("googl_token", "bitly_token"))

# Bitly_api_version <- "v3"
# Googl_api_version <- "v1"
# Isgd_api_version <- "v2015"

#' @title Assign API tokens using OAuth2.0
#' 
#' @description You should register an application in order to get Client ID and Client Secret code. 
#' For Bit.ly, go to \url{https://bitly.com/a/oauth_apps} and in the field \code{Redirect URIs:} 
#' type for example "http://localhost:1410". 
#' For Goo.gl API Keys you should go to the \code{http://console.developers.google.com/project/},
#' select "APIs & auth", then "Credentials", then "add OAuth2.0 client ID" and lastly you select 
#' "Type:Other". 
#' 
#' @param key - Client ID
#' @param secret - Client Secret
#' 
#' @seealso See \url{http://dev.bitly.com/rate_limiting.html}
#' @seealso See \url{http://dev.bitly.com/authentication.html}
#' @seealso See \url{https://developers.google.com/url-shortener/v1/getting_started#APIKey}
#' 
#' @examples
#' \dontrun{
#' googl_token <-
#'   googl_auth(key = "806673580943-78jdskus76fu7r0m21erihqtltcka29i.apps.googleusercontent.com",
#'              secret = "qItL-PZnm8GFxUOYM0zPVr_t")
#' bitly_token <-
#'   bitly_auth(key = "be03aead58f23bc1aee6e1d7b7a1d99d62f0ede8",
#'              secret = "b7e4abaf8b26ec4daa92b1e64502736f5cd78899")
#' }
#' 
#' @import httr
#' @export
googl_auth <- function(key = "", secret = "") {
  googl_token <- httr::oauth2.0_token(httr::oauth_endpoints("google"),
                                       httr::oauth_app("google", key = key, secret = secret),
                                       scope = "https://www.googleapis.com/auth/urlshortener",
                                       cache = TRUE)
  .state$token <- googl_token
  return(googl_token)
}


#' @rdname googl_auth
#' @export
bitly_auth <- function(key = "", secret = "") {
  bitly_token <- httr::oauth2.0_token(httr::oauth_endpoint(authorize = "https://bitly.com/oauth/authorize",
                                                           access = "https://api-ssl.bitly.com/oauth/access_token"),
                                      httr::oauth_app("bitly", key = key, secret = secret),
                                      cache = TRUE)
  .state$token <- bitly_token
  return(bitly_token)
}

#' @title Generalized function for executing GET/POST requests
#' 
#' @param url - which is used for the request
#' @param authcode - calls the rbitlyApi \code{\link{rbitlyApi}}
#' @param queryParameters - parameters that are used for building a URL
#' @param showURL - for debugging purposes only: it shows what URL has been called
#' 
#' @import httr
#' @import jsonlite
#' 
#' @noRd
#' @keywords internal
doRequest <- function(verb, url, service = "", queryParameters = NULL, showURL = NULL) {
  service_token <- if (service == "bitly") {
    bitly_token
  }
  
  service_token <- if (service == "googl") {
    googl_token
  } else {
    NULL
  }
  
  switch(verb,
         "GET" = {
           return_request <- httr::GET(url, query = queryParameters, httr::config(token = service_token))
         },
         "POST" = {
           return_request <- httr::POST(url, body = queryParameters, encode = "json", 
                                        httr::content_type_json(), httr::config(token = service_token))
         }
  )
  
  if (http_error(return_request) == FALSE) {
    text_response <- content(return_request, as = "text")
    json_response <- fromJSON(text_response)
    
    if (is.null(json_response$status_code) == FALSE && json_response$status_code >= 400) {
      message(sprintf("Code: %s - %s", json_response$status_code, json_response$status_txt))
    }
      
    if (identical(showURL, TRUE)) {
      cat("The requested URL has been this: ", return_request$request$url, "\n") 
    }
    
  } else {
    stop_for_status(return_request)
  }
  
  return(json_response)
}



