##' @include methods.R
NULL


##' Basic timer widget
##'
##' Calls FUN every ms/1000 seconds. A timer is stopped through its \code{stop_timer} method which is called using OO style: \code{obj$stop_timer()}.
##' @param ms interval in milliseconds
##' @param FUN FUnction to call. Has one argument, data passed in
##' @param data passed to function
##' @param one.shot logical. If TRUE, called just once, else repeats
##' @param start logical. If FALSE, started by \code{start_timer} OO method. (Call \code{obj$start_time()}). 
##' @param toolkit gui toolkit to dispatch into
##' @export
##' @rdname gtimer
##' @examples
##' \dontrun{
##' i <- 0
##' FUN <- function(data) {i <<- i + 1; if(i > 10) a$stop_timer(); print(i)}
##' a <- gtimer(1000, FUN)
##' ##
##' ## a one shot timer is run only once
##' FUN <- function(data) message("Okay, I can breathe now")
##' hold_breath <- gtimer(1000*60, FUN, one.shot=TRUE)
##' }
gtimer <- function(ms, FUN, data=NULL,  one.shot=FALSE, start=TRUE, toolkit=guiToolkit()) {
  .gtimer(toolkit=toolkit, ms=ms, FUN=FUN, data=data, one.shot=one.shot, start=start)
}

##' S3 generic for dispatch
##'
##' @export
##' @rdname gtimer
.gtimer <- function(toolkit, ms, FUN, data=NULL, one.shot=FALSE, start=TRUE) UseMethod(".gtimer")
