##' @include methods.R
NULL

##' Basic progress bar widget
##'
##' @inheritParams gwidget
##' @return a \code{GButton} instance. While this object has its own
##' (reference) methods, one primarily interacts with it through S3
##' methods defined within the package.
##' @param value Initial value, between 0 and 100. A value of \code{NULL} will make pulsing bar with indeterminate state. For some toolkits, this must be called periodically to pulse the bar.
##' @export
##' @examples
##' \dontrun{
##' w <- gwindow("progress bar example")
##' pb <- gprogressbar(cont=w)
##' for(i in 10:100) {Sys.sleep(.1); svalue(pb) <- i}
##' }
gprogressbar <- function(value=10, container=NULL, ..., toolkit=guiToolkit()) {

  if(is.character(toolkit))
    toolkit <- guiToolkit(toolkit)

  
  obj <- .gprogressbar(toolkit, value, container, ...)

  check_return_class(obj, "GProgressBar")
  obj
  
}

##' S3 generic whose methods are implemented in the toolkit packages
##'
##' @rdname gprogressbar
##' @export
.gprogressbar <- function(toolkit, value, container, ...) UseMethod(".gprogressbar")

