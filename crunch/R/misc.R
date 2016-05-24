is.error <- function (x) inherits(x, "try-error")

rethrow <- function (x) halt(errorMessage(x))

errorMessage <- function (e) attr(e, "condition")$message

updateList <- function (x, y) {
    x[names(y)] <- y
    return(x)
}

##' Generic List Element Extractor
##'
##' @param key character naming the key(s) to extract. Can traverse list
##' elements by separating them with \code{$}.
##' @param xlist list containing other lists from which you want to extract
##' @param ifnot what to return if the key is not found in a given xlist element
##' @param simplify logical, passed to sapply internally
##' @return the requested element(s). If length(key)>1, a named list of those
##' elements
##' @keywords internal
selectFrom <- function (key, xlist, ifnot=NA, simplify=TRUE) {
    if (!is.list(xlist)) {
        halt("xlist must be a list object")
    }
    if (length(key)>1) {
        y <- sapply(key, selectFrom, xlist, ifnot, simplify=FALSE)
    } else {
    	y <- sapply(xlist,
    	    function (x) {
    	        key <- unlist(strsplit(key, "$", fixed=TRUE))
    	        for (i in key) {
    	            if (!is.list(x)) x <- NULL
                    if (!is.null(x)) x <- x[[i]]
                }
                if (is.null(x)) x <- ifnot
                return(x)
    	    }, simplify=simplify)
    }
    return(y)
}

vget <- function (name) {
    ## Return a function you can lapply/vapply to select an attribute
    ## Usage: lapply(list.of.stuff, vget("name"))
    ## instead of: lapply(list.of.stuff, function (x) x$name)
    return(function (x) x[[name]])
}

##' Make a prose list
##'
##' Function to paste together a list of items, separated by commas (if more
##' than 2), and with the last one having the collapse string.
##'
##' @param x vector or list
##' @param collapse default="and"
##' @keywords internal
serialPaste <- function (x, collapse="and") {
	if (length(x)>1) x[length(x)] <- paste(collapse, x[length(x)])
	join.with <- ifelse(length(x)>2, ", ", " ")
	return(paste(x, collapse=join.with))
}

now <- function () strftime(Sys.time(), usetz=TRUE)

##' @importFrom httr parse_url build_url
absoluteURL <- function (urls, base) {
    ## Detect if we have relative urls, and then concatenate if so
    if (length(urls) && ## if there is anything to munge
        !any(substr(urls, 1, 4) == "http")) { ## the urls don't start with http
            base.url <- parse_url(base)
            urls <- vapply(urls, function (x, b) {
                b$path <- joinPath(b$path, x)
                if (is.null(b$scheme)) {
                    ## If file path and not URL, as in for tests,
                    ## let's return it relative
                    return(b$path)
                }
                ## Pop off any leading "/" because build_url will add it
                b$path <- sub("^/", "", b$path)
                b$query <- NULL ## Catalog query params aren't valid for entities
                return(build_url(b))
            }, character(1), b=base.url, USE.NAMES=FALSE)
        }
    return(urls)
}

joinPath <- function (base.path, relative.part) {
    first.char <- substr(relative.part, 1, 1)
    if (first.char == "/") {
        ## This is absolute, relative to the host
        return(relative.part)
    }
    u <- c(strsplit(base.path, "/")[[1]], strsplit(relative.part, "/")[[1]])
    ## Drop any references to current location (.)
    u <- u[u != "."]
    ## Walk the ..
    if (any(u == "..")) {
        ## If we're here, we must have some normalization to do
        i <- 1
        n <- length(u)
        while (i <= n) {
            if (u[i] == "..") {
                ## Remove i and the one before it, and roll the counter back
                u <- u[-c(i-1, i)]
                n <- n - 2
                i <- i - 1
            } else {
                i <- i + 1
            }
        }
    }
    out <- paste(u, collapse="/")
    last.char <- substr(relative.part, nchar(relative.part),
        nchar(relative.part))
    if (last.char == "/") {
        out <- paste0(out, "/")
    }
    return(out)
}

askForPermission <- function (prompt="") {
    ## If options explicitly say we don't need to ask, bail.
    ## Have to check that it's FALSE and not NULL. Silence doesn't mean consent.
    must.confirm <- getOption("crunch.require.confirmation") %||% TRUE
    if (must.confirm == FALSE) return(TRUE)

    ## If we're here but not interactive, we can't give permission.
    if (!interactive()) return(FALSE)
    prompt <- paste(prompt, "(y/n) ")
    proceed <- ""
    while (!(proceed %in% c("y", "n"))) {
        proceed <- tolower(readline(prompt))
    }
    return(proceed == "y")
}

emptyObject <- function () {
    ## toJSON(list()) is "[]". toJSON(emptyObject()) is "{}"
    structure(list(), .Names=character(0))
}

## Borrowed from Hadley
"%||%" <- function (a, b) if (!is.null(a)) a else b
