##' @include methods.R
NULL

##'  Spinbutton constructor
##'
##' A spinbutton allows the user to select from a pre-selected range
##' of numbers. Similar to a slider, but with more precision, but
##' slower to adjust. The basic arguments mirror that of \code{seq.int}.
##' @param from from value
##' @param to to value
##' @param by step length
##' @param length.out number of steps. Only one of \code{by} or \code{length.out} is used.
##' @param along.with Take from
##' @param value initial value
##' @param digits number of digits to display, should the toolkit support it
##' @inheritParams gwidget
##' @export
##' @rdname gspinbutton
##' @seealso \code{\link{gslider}}
##' @example inst/examples/ex-rangewidget.R
gspinbutton =function(
  from = 0, to = 10, by = 1,
  length.out = NULL, along.with=NULL,
  value = from, digits = 0,
  handler = NULL, action = NULL, container = NULL, ... ,
  toolkit=guiToolkit()){

  ## mostly from seq.default
  if (!missing(along.with)) {
    length.out <- length(along.with)
  } else if(!missing(length.out)) {
    len <- length(length.out)
    if (!len) 
      stop("argument 'length.out' must be of length 1")
    if (len > 1L) {
      warning("first element used of 'length.out' argument")
      length.out <- length.out[1L]
    }
    length.out <- ceiling(length.out)
  }
  if(!is.null(length.out)) {
    ## set up by to be for length,out
    by <- (to - from)/(length.out[1] - 1)
  }

  obj <- .gspinbutton(toolkit, from, to, by, value, digits,
                      handler, action, container=container, ...)
  check_return_class(obj, "GSpinButton")
  obj
  
}

##' generic for toolkit dispatch
##'
##' @export
##' @rdname gspinbutton
.gspinbutton <- function(toolkit,
                         from = 0, to = 10, by = 1, value = from, digits = 0,
                         handler = NULL, action = NULL, container = NULL, ... ) UseMethod( '.gspinbutton' )
