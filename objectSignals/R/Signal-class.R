## ======================================================================
## Signal class
## ======================================================================
Signal.gen <- setRefClass("Signal",
                          fields = list(.listeners = "list", emit = "function",
                            .idCounter = "integer", .blocked = "logical",
                            .buffered = "logical", .queue = "list",
                            .accumulator = "function"))

##' Creates a \code{Signal} object, with which a mutable object can report
##' changes to its state. Interested clients register a function that
##' is called whenever the signal is emitted. This allows those
##' clients to respond to changes in the object state.
##'
##' A \code{Signal} object is usually created by a constructor and
##' stored as a field in a reference class object. Clients then access
##' the signal either directly or through an accessor.
##'
##' The \code{Signal} reference class has the following methods:
##' \describe{
##'   \item{connect(FUN, namedArgs = FALSE, ...)}{Connect \code{FUN} as a
##'     handler for the signal. A unique identifier is returned, which
##'     can be used to later disconnect the handler.
##'     Handler invocation follows these rules:
##'     \itemize{
##'       \item{\code{namedArgs=TRUE} arguments are named in call to handler.
##'         Otherwise, they are unnamed and matching is by position.}
##'       \item{If a handler is missing a signal argument, the argument is
##'         dropped when calling the handler.}
##'       \item{A handler may have arguments not in the signal signature.}
##'       \item{Arguments in \code{...} are appended to the handler call.}
##'     }
##'   }
##'   \item{disconnect(id)}{Disconnects the handler registered with the
##'     identifier \code{id}.}
##'   \item{emit(<signal signature>)}{Emits the signal, calling all of its
##'     handlers with the passed arguments. The signature depends on how the
##'     signal was constructed. All signal args must be passed to \code{emit},
##'     unless they have a default.
##'   }
##'   \item{block()}{Blocks emission of the signal. All emission requests are
##'     ignored.}
##'   \item{unblock()}{Unblock the signal.}
##'   \item{buffer()}{Buffer emissions, waiting to pass them to the handlers
##'    until \code{flush} is called.}
##'   \item{flush()}{Flush the emission buffer, calling every handler on every
##'    buffered emission.}
##'   \item{accumulator(value)}{If called with no arguments, get the function,
##'    if any, used to combine events in the buffer into a single event.
##'    Otherwise, \code{value} replaces the current function. The accumulator
##'    function should take one or two arguments. If it takes one argument, it
##'    is invoked upon a flush and is passed the list of events in the buffer.
##'    An event is simply a list containing the arguments passed to \code{emit}.
##'    If the accumulator function takes two arguments, the function is invoked
##'    upon every emission, when buffering is active and there is one event in
##'    the buffer. The first argument is the currently buffered event and the
##'    second is the new event that the function should merge into the first.
##'    The returned event then replaces the event in the buffer.}   
##' }
##' @title Signals
##' @param ... Arguments that express the signature of the signal.
##' @return An instance of the reference class \code{Signal}
##' @author Michael Lawrence
##' @aliases Signal-class
##' @example objectSignals/inst/examples/Signal.R
Signal <- function(...) {
  call <- sys.call()[-1L]
  hasDefault <-
    if (is.null(names(call)))
      rep(FALSE, length(call))
    else nzchar(names(call))
  names(call)[!hasDefault] <- sapply(call[!hasDefault], deparse)
  call[!hasDefault] <- alist(foo=)
  signal <- Signal.gen$new(.idCounter = 0L, .blocked = FALSE, .buffered = FALSE)
  signal$emit <- as.function(c(as.list(call), quote({
    if (.blocked)
      return(invisible())
    if (.buffered) {
      event <- lapply(match.call()[-1], eval, environment())
      if (length(formals(.accumulator)) == 2L && length(.queue))
        .queue[[1]] <<- .accumulator(.queue[[1]], event)
      else .queue <<- c(.queue, list(event))
    } else {
      for (listener in .listeners)
        eval(listener)
    }
    invisible()
  })), signal@.xData) # is getting .xData bad practice?
  signal
}

Signal.gen$methods(connect = function(FUN, namedArgs = FALSE, ...) {
  ## FUN is a function
  extraArgs <- list(...)
  if (length(extraArgs) > length(formals(FUN)))
    stop("More args in '...' than in formals of 'FUN'")
  args <- sapply(names(formals(emit)), as.name)
  ## the wrapper adds the extra args, and also uses '...' to drop unwanted args
  wrapperFormals <- formals(FUN)
  handlerArgs <- sapply(names(wrapperFormals), as.name)
  if (!("..." %in% names(wrapperFormals)))
    wrapperFormals <- c(wrapperFormals, alist(...=))
  wrapper <- as.function(c(wrapperFormals, as.call(c(FUN, handlerArgs))))
  .idCounter <<- .idCounter + 1L
  id <- as.character(.idCounter)
  if (!namedArgs)
    names(args) <- NULL
  .listeners[[id]] <<- as.call(c(list(wrapper), c(args, extraArgs)))
  invisible(id)
})

Signal.gen$methods(disconnect = function(id) {
  .listeners[[id]] <<- NULL
  invisible(.self)
})

Signal.gen$methods(block = function() {
  .blocked <<- TRUE
  invisible(.self)
})

Signal.gen$methods(unblock = function() {
  .blocked <<- FALSE
  invisible(.self)
})

Signal.gen$methods(buffer = function() {
  .buffered <<- TRUE
  invisible(.self)
})

Signal.gen$methods(flush = function() {
  if (length(formals(.accumulator)) == 1L)
    .queue <<- list(.accumulator(.queue))
  .buffered <<- FALSE
  .accumulator <<- function() NULL
  lapply(.queue, do.call, what = emit)
  .queue <<- list()
  invisible(.self)
})

## Allows C++-style initializer chaining
Signal.gen$methods(accumulator = function(value) {
  if (missing(value))
    .accumulator
  else {
    .accumulator <<- value
    invisible(.self)
  }
})

setMethod("show", "Signal", function(object) {
  cat(deparse(as.call(c(as.name(class(object)), formals(object$emit)))), "with",
      length(object$.listeners), "listeners\n")
})

##' Return how many listerns a \code{Signal} object contain
##'
##' @title length method for \code{Signal} class
##' @param x \code{Signal} object
##' @return Numieric value indicate how many listeners this signal contains
##' @rdname length-methods
##' @seealso \code{\link{listeners}}
##' @author Tengfei Yin
setMethod("length", "Signal", function(x){
  length(x$.listeners)
})

setGeneric("listeners", function(object) standardGeneric("listeners"))

##' Return listerners of a \code{Signal} object
##'
##' @title listeners method for \code{Signal} class
##' @param object \code{Signal} object
##' @return A list of listeners
##' @docType methods
##' @rdname listeners-methods
##' @author Tengfei Yin 
setMethod("listeners", "Signal", function(object){
  object$.listeners
})


