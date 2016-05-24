
## TODO: possibly allow a 'slot' parameter, and only pass the
## indicated arguments. 
qconnect <- function(x, signal, handler, user.data)
{
  stopifnot(is(x, "QObject"))
  stopifnot(is.function(handler))
  signal <- as.character(signal)
  has.user.data <- !missing(user.data)
  if (!has.user.data)
    user.data <- NULL
  nargs <- length(formals(handler))
  if (has.user.data && !nargs)
    stop("if 'user.data' specified, handler must take at least one argument")
  hasDots <- "..." %in% names(formals(handler))
  if (hasDots && nargs == 1) # only '...'
    ## could be a wrapper, so we do not know actual argument count
    signal <- qresolveSignature(x, signal, "signal")
  else {
    if (hasDots)
      nargs <- Inf
    signal <- qresolveSignature(x, signal, "signal", nargs)
  }
  .Call("qt_qconnect", x, signal, handler, user.data, has.user.data, PACKAGE="qtbase")
}
