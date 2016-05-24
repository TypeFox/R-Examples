#' R interface for Google Civic API
#'
#' See \url{(https://developers.google.com/civic-information/)}
#'
#' @import jsonlite RCurl
#' @export
#' @param method method to call
#' @param key API key
#' @param debug Enable debugging mode
#' @keywords API
#' @seealso \code{\link{PopongAPI}}, \code{\link{SunlightAPI}}
#'

GoogleAPI <- function(method, key=getOption("GoogleAPIKey"), debug=FALSE) {
    # TODO: auto navigate pages
    apiSource   <- "google"
    apiAttrs    <- eval(parse(text=sprintf("apiInfo$%s", apiSource)))
    apiVersion  <- paste0("v", as.integer(apiAttrs$version))

    MethodInAPI(apiSource, method)

    paths <- c(apiAttrs$url, apiVersion, method)
    query  <- list("key"=key)
    url   <- FormatURL(paths, query)
    jsontext <- RCurl::getURL(url)
    response <- jsonlite::fromJSON(jsontext)
    if(debug) { print(url) }

    return(response)
}
