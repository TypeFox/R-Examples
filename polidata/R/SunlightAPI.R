#' R interface for Sunlight API
#'
#' See \url{(https://sunlightlabs.github.io/congress/)}
#'
#' @import jsonlite
#' @export
#' @param method method to call
#' @param key API key
#' @param debug Enable debugging mode
#' @keywords API
#' @seealso \code{\link{GoogleAPI}}, \code{\link{PopongAPI}}
#'

SunlightAPI <- function(method, key=getOption("SunlightAPIKey"), debug=FALSE) {
    # TODO: auto navigate pages
    apiSource   <- "sunlight"
    apiAttrs    <- eval(parse(text=sprintf("apiInfo$%s", apiSource)))

    MethodInAPI(apiSource, method)

    paths <- c(apiAttrs$url, method)
    query <- list("apikey"=key)
    url   <- FormatURL(paths, query)
    jsontext <- RCurl::getURL(url)
    response <- jsonlite::fromJSON(jsontext)
    if(debug) { print(url) }

    return(response$results)
}
