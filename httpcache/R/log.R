##' Log a message
##'
##' @param ... Strings to pass to \code{cat}
##' @return Nothing
##' @export
logMessage <- function (...) {
    logfile <- getOption("httpcache.log")
    if (!is.null(logfile)) {
        cat(strftime(Sys.time(), "%Y-%m-%dT%H:%M:%S"), ..., "\n", file=logfile, append=TRUE)
    }
}

##' Stop, log, and no call
##'
##' Wrapper around \code{stop} that logs the error message and then stops
##' with call.=FALSE by default.
##' @param ... arguments passed to \code{stop}
##' @param call. logical: print the call? Default is \code{FALSE}, unlike
##' \code{stop}
##' @return Nothing. Raises an error.
##' @export
halt <- function (..., call.=FALSE) {
    msg <- gsub("\n", " ", ..1)
    logMessage("ERROR", msg)
    stop(..., call.=call.)
}

##' Enable logging
##'
##' @param filename character: a filename/path where the log can be written out.
##' If \code{""}, messages will print to stdout (the screen). See
##' \code{\link[base]{cat}}.
##' @param append logical: if the file already exists, append to it? Default
##' is \code{FALSE}, and if not in append mode, if the \code{filename} exists,
##' it will be deleted.
##' @return Nothing.
##' @export
startLog <- function (filename="", append=FALSE) {
    options(httpcache.log=filename)
    if (!append && nchar(filename) && file.exists(filename)) {
        file.remove(filename)
    }
}

##' Read in a httpcache log file
##'
##' @param filename character name of the log file, passed to
##' \code{\link[utils]{read.delim}}
##' @param scope character optional means of selecting only certain log
##' messages. By default, only "CACHE" and "HTTP" log messages are kept. Other
##' logged messages, such as "ERROR" messages from \code{\link{halt}}, will be
##' dropped from the resulting data.frame.
##' @return A data.frame of log results.
##' @export
##' @importFrom utils read.delim
loadLogfile <- function (filename, scope=c("CACHE", "HTTP")) {
    df <- read.delim(filename, sep=" ", header=FALSE,
        stringsAsFactors=FALSE)[,1:6]
    names(df) <- c("timestamp", "scope", "verb", "url", "status", "time")
    df <- df[df$scope %in% scope,] ## Prune out-of-scope things
    df$timestamp <- strptime(df$timestamp, "%Y-%m-%dT%H:%M:%S")
    df[c("status", "time")] <- lapply(df[c("status", "time")], as.numeric)
    return(df)
}

##' Summarize cache performance from a log
##'
##' @param logdf A logging data.frame, as loaded by \code{link{loadLogfile}}.
##' @return A list containing counts of cache hit/set/drop events, plus a
##' cache hit rate.
##' @export
cacheLogSummary <- function (logdf) {
    df <- logdf[logdf$scope == "CACHE",]
    counts <- table(df$verb)
    return(list(counts=counts,
        hit.rate=100*counts["HIT"]/sum(counts[c("HIT", "SET")])))
}

##' Summarize HTTP requests from a log
##'
##' @param logdf A logging data.frame, as loaded by \code{link{loadLogfile}}.
##' @return A list containing counts of HTTP requests by verb, as well as
##' summaries of time spent waiting on HTTP requests.
##' @export
##' @importFrom utils head tail
requestLogSummary <- function (logdf) {
    total.time <- as.numeric(difftime(tail(logdf$timestamp, 1),
        head(logdf$timestamp, 1), units="secs"))
    df <- logdf[logdf$scope == "HTTP",]
    counts <- table(df$verb)
    req.time <- sum(df$time, na.rm=TRUE)
    pct.http.time <- 100*req.time/total.time
    return(list(counts=counts, req.time=req.time, total.time=total.time,
        pct.http.time=pct.http.time))
}
