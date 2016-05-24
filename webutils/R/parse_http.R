#' Parse http request
#'
#' Parse the body of a http request, based on the \code{Content-Type} request
#' header. Currently supports the three most important content types:
#' \code{application/x-www-form-urlencoded} (\code{\link{parse_query}}),
#' \code{multipart/form-data} (\code{\link{parse_multipart}}) and
#' \code{application/json} (\code{\link{fromJSON}}).
#'
#' @export
#' @param body request body of the http request
#' @param content_type content-type http request header as specified by the client
#' @param ... additional arguments passed to parser function
#' @importFrom jsonlite fromJSON
#' @examples # Parse json encoded payload:
#' parse_http('{"foo":123, "bar":true}', 'application/json')
#'
#' # Parse url-encoded payload
#' parse_http("foo=1%2B1%3D2&bar=yin%26yang", "application/x-www-form-urlencoded")
#'
#' \dontrun{use demo app to parse multipart/form-data payload
#' demo_rhttpd()
#' }
parse_http <- function(body, content_type, ...){
  # Remove header name if present
  content_type <- sub("Content-Type: ?", "", content_type, ignore.case=TRUE);

  # Switch by content-type
  if(grepl("multipart/form-data; boundary=", content_type, fixed=TRUE)){
    return(parse_multipart(body, get_boundary(content_type)))
  } else if(grepl("application/x-www-form-urlencoded", content_type, fixed=TRUE)){
    return(parse_query(body))
  } else if(grepl("(text|application)/json", content_type)){
    if(is.raw(body))
      body <- rawToChar(body)
    return(fromJSON(body, ...))
  } else {
    stop("Unsupported Content-Type: ", content_type)
  }
}
