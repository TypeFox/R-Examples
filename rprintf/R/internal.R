setnames <- `names<-`

patterns <- c(
  rsprintf = "(?<!%)%[\\s\\+\\-\\#\\.\\d]*\\w",
  rprintv = "(?<!\\$)\\$[\\w\\.]+(:[\\s\\+\\-\\#\\.\\d]*\\w)?",
  rprintn = "(?<!\\{)\\{\\d+(:[\\s\\+\\-\\#\\.\\d]*\\w)?\\}(?!\\})",
  rsprintf = ".")

rprintf.match <- function(x, fun, args) {
  do.call(fun, c(x, args))
}

makelist <- function(...) {
  ## better use do.call(c,list(...))
  args <- list(...)
  if (any(vapply(args, is.list, logical(1L))))
    unlist(args, recursive = FALSE) else args
}

rsprintf <- function(x, ...) {
  args <- makelist(...)
  do.call(sprintf, c(x, args))
}
