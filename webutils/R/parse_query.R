#' Parse query string
#'
#' Parse http parameters from a query string. This includes unescaping
#' of url-encoded values.
#'
#' For http GET requests, the query string is specified
#' in the URL after the question mark. For http POST or PUT requests, the query
#' string can be used in the request body when the \code{Content-Type} header
#' is set to \code{application/x-www-form-urlencoded}.
#'
#' @export
#' @importFrom utils URLdecode
#' @param query a url-encoded query string
#' @examples q <- "foo=1%2B1%3D2&bar=yin%26yang"
#' parse_query(q)
parse_query <- function(query){
  if(is.raw(query)){
    query <- rawToChar(query);
  }
  stopifnot(is.character(query));

  #httpuv includes the question mark in query string
  query <- sub("^[?]", "", query)

  #split by & character
  argslist <- sub("^&", "", regmatches(query, gregexpr("(^|&)[^=]+=[^&]+", query))[[1]])
  argslist <- strsplit(argslist, "=");
  values <- lapply(argslist, function(x){if(length(x) < 2) "" else paste(x[-1], collapse="=")});
  values <- lapply(values, function(s) {URLdecode(chartr('+',' ',s))});
  names(values) <- lapply(argslist, "[[", 1);
  return(values)
}
