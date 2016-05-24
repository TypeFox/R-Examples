##' Cache-aware versions of httr verbs
##'
##' These functions set, read from, and bust the HTTP query cache. They wrap
##' the similarly named functions in the httr package and can be used as
##' drop-in replacements for them.
##'
##' \code{GET}
##' checks the cache before making an HTTP request, and if there is a cache
##' miss, it sets the response from the request into the cache for future
##' requests. The other verbs, assuming a more or less RESTful API, would be
##' assumed to modify server state, and thus they should trigger cache
##' invalidation. They have default cache-invalidation strategies, but you can
##' override them as desired.
##'
##' @param url character URL of the request
##' @param ... additional arguments passed to the httr functions
##' @param drop For \code{PUT}, \code{PATCH}, \code{POST}, and \code{DELETE},
##' code to be executed after the request. This is intended to be for supplying
##' cache-invalidation logic. By default, \code{POST} drops cache only for
##' the specified \code{url} (i.e. \code{\link{dropOnly}}), while the other
##' verbs drop cache for the request URL and for any URLs nested below it
##' (i.e. \code{\link{dropCache}}).
##' @return The corresponding httr response object, potentially read from cache
##' @importFrom httr GET
##' @importFrom digest digest
##' @aliases GET PUT POST PATCH DELETE
##' @name cached-http-verbs
##' @seealso \code{\link{dropCache}}
##' @export
GET <- function (url, ...) {
    if (caching()) {
        Call <- match.call(expand.dots = TRUE)
        cache.url <- url
        if (!is.null(Call[["query"]])) {
            cache.url <- paste0(url, "?HASHED_QUERY=",
                digest(eval.parent(Call$query)))
        }
        if (exists(cache.url, envir=cache)) {
            logMessage("CACHE HIT", cache.url)
            return(get(cache.url, envir=cache))
        }
    }
    x <- httr::GET(url, ...)
    logMessage(responseStatusLog(x))
    if (caching() && x$status_code == 200) {
        logMessage("CACHE SET", cache.url)
        assign(cache.url, x, envir=cache)
    }
    return(x)
}

##' @rdname cached-http-verbs
##' @export
##' @importFrom httr PUT
PUT <- function (url, ..., drop=dropCache(url)) {
    x <- httr::PUT(url, ...)
    logMessage(responseStatusLog(x))
    force(drop)
    return(x)
}

##' @rdname cached-http-verbs
##' @export
##' @importFrom httr POST
POST <- function (url, ..., drop=dropOnly(url)) {
    x <- httr::POST(url, ...)
    logMessage(responseStatusLog(x))
    force(drop)
    return(x)
}

##' @rdname cached-http-verbs
##' @export
##' @importFrom httr PATCH
PATCH <- function (url, ..., drop=dropCache(url)) {
    x <- httr::PATCH(url, ...)
    logMessage(responseStatusLog(x))
    force(drop)
    return(x)
}

##' @rdname cached-http-verbs
##' @export
##' @importFrom httr DELETE
DELETE <- function (url, ..., drop=dropCache(url)) {
    x <- httr::DELETE(url, ...)
    logMessage(responseStatusLog(x))
    force(drop)
    return(x)
}

responseStatusLog <- function (response) {
    req <- response$request
    return(paste("HTTP",
        req$method,
        req$url,
        response$status_code,
        response$times["total"]))
}
