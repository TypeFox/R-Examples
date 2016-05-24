#' This will request (POST method) with the prepared Query to the Google Analytics
#' and returns the response in R object (res).
#'
#' @keywords internal
#'
#' @param query.uri The URL to insert report in API.
#' @param postbody  Body part of POST request.
#' @param token Oauth2.0 token.
#' @return res Response of POST request
#'
#' @importFrom httr POST
#' @importFrom httr content_type_json
#' @importFrom httr add_headers
#'
GetPostResponse <- function(query.uri,postbody,token) {

  tok <- paste("Bearer", token$credentials$access_token)
  res <- POST(url = query.uri, body = postbody, content_type_json(), add_headers(Authorization = tok))

  return(res)
}
