base64encode <- function(what, linewidth, newline) {
  linewidth <- if (missing(linewidth) || !is.numeric(linewidth) || length(linewidth) < 1L) 0L else as.integer(linewidth[1L])
  if (is.na(linewidth)) linewidth <- 0L else if (linewidth > 0L && linewidth < 4L) linewidth <- 4L
  if (missing(newline)) newline <- NULL  
  fi <- NULL
  if (is.character(what)) {
    what <- file(what, "rb")
    on.exit(close(what))
  }
  if (inherits(what, "connection")) {
    slice <- 65535L  ## default slice size - must be divisible by 3
    if (linewidth > 0L) { ## we have to make sure the slices span whole lines
      if (linewidth %% 4L > 0) linewidth <- linewidth - linewidth %% 4L
      bw <- as.integer(linewidth / 4L) * 3L
      if (slice %% bw > 0L)
        slice <- slice + (bw - (slice %% bw))
    }
    l <- list()
    while (length(r <- readBin(what, raw(0), slice)))
      l <- c(l, .Call(B64_encode, r, linewidth, newline))
    if (linewidth > 0L && is.null(newline))
      unlist(l)
    else paste(unlist(l), collapse = if (is.null(newline)) "" else newline)
  } else
  .Call(B64_encode, as.raw(what), linewidth, newline)
}

base64decode <- function(what, output=NULL, file) {
  if (!missing(file) && !missing(what)) stop("'what' and 'file' are mutually exclusive")
  if (!missing(file)) {
    what <- file(file, "r")
    on.exit(close(what))
  }
  if (is.character(output)) {
    output <- file(output, "wb")
    on.exit(close(output))
  } else if (!inherits(output, "connection") && !is.null(output)) stop("output must be a filename, connection or NULL")
  r <- if (inherits(what, "connection")) {
    ## FIXME: we may want to use chunking ...
    .Call(B64_decode, readLines(what, warn=FALSE))
  } else
    .Call(B64_decode, what)
  if (inherits(output, "connection")) {
    writeBin(r, output)
    invisible(length(r))
  } else r
}
