## Create the cache env
cache <- NULL
initCache <- function () {
    cache <<- new.env(hash=TRUE)
}
initCache()

caching <- function () {
    ## Default should be on, so if httpcache.on isn't set, return TRUE
    opt <- getOption("httpcache.on")
    return(is.null(opt) || isTRUE(opt))
}
##' Manage the HTTP cache
##'
##' These functions turn the cache on and off and clear the contents of the
##' query cache.
##' @return Nothing. Functions are run for their side effects.
##' @aliases cacheOn cacheOff clearCache
##' @name cache-management
##' @export
cacheOn <- function () options(httpcache.on=TRUE)

##' @rdname cache-management
##' @export
cacheOff <- function () {
    options(httpcache.on=FALSE)
    clearCache()
}

##' @rdname cache-management
##' @export
clearCache <- function () {
    logMessage("CACHE CLEAR")
    rm(list=ls(all.names=TRUE, envir=cache), envir=cache)
}

##' Context manager to temporarily turn cache off if it is on
##'
##' If you don't want to store the response of a GET request in the cache,
##' wrap it in \code{uncached()}. It will neither read from nor write to cache.
##' However, \code{uncached} will not invalidate cache records, if present.
##'
##' @param ... Things to evaluate with caching off
##' @return Whatever ... returns.
##' @examples
##' uncached(GET("http://httpbin.org/get"))
##' @export
uncached <- function (...) {
    old <- getOption("httpcache.on")
    on.exit(options(httpcache.on=old))
    options(httpcache.on=FALSE)
    eval.parent(...)
}

##' Invalidate cache
##'
##' These functions let you control cache invalidation. \code{dropOnly}
##' invalidates cache only for the specified URL. \code{dropPattern} uses
##' regular expression matching to invalidate cache. \code{dropCache} is a
##' convenience wrapper around \code{dropPattern} that invalidates cache for
##' any resources that start with the given URL.
##' @param x character URL or regular expression
##' @return Nothing. Functions are run for their side effects.
##' @export
dropCache <- function (x) {
    ## Drop x and anything below it in the tree
    dropPattern(paste0("^", regexEscape(popQuery(x))))
}

##' @rdname dropCache
##' @export
dropOnly <- function (x) {
    logMessage("CACHE DROP", x)
    suppressWarnings(rm(list=x, envir=cache))
}

##' @rdname dropCache
##' @export
dropPattern <- function (x) {
    logMessage("CACHE DROP", x)
    rm(list=ls(envir=cache, pattern=x), envir=cache)
}

# dropBelow <- function (x) {
#     ## Don't drop x, just those below it in the tree. hence ".+"
#     dropPattern(paste0("^", regexEscape(popQuery(x)), ".+"))
# }

regexEscape <- function (x) {
    ## Escape all reserved characters that are valid URL chars with \\
    for (i in unlist(strsplit(".+?*", ""))) {
        x <- gsub(paste0("(\\", i, ")"), "[\\1]", x)
    }
    return(x)
}

popQuery <- function (x) {
    ## Remove query parameters from a URL
    return(sub("\\?.*$", "", x))
}

.internalFunction <- function () TRUE
