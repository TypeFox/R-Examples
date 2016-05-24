.log.info <- new.env()

str2int <- function(str)
  switch(str, error=1, warn=2, info=3, debug=4, verbose=5, 5)

loglevel <- function(lev) {
  lev <- deparse(substitute(lev))
  .log.info$level <- str2int(lev)
  invisible(.log.info$level)
}

verbose <- function(...) {
  if (!is.null(.log.info$level) && .log.info$level >= str2int('verbose')) {
    msg <- paste('verbose debug: ', paste(list(...), collapse=''), '\n', sep='')
    cat(msg)
  }
}

dbug <- function(...) {
  if (!is.null(.log.info$level) && .log.info$level >= str2int('debug')) {
    msg <- paste('debug: ', paste(list(...), collapse=''), '\n', sep='')
    cat(msg)
  }
}

info <- function(...) {
  if (!is.null(.log.info$level) && .log.info$level >= str2int('info')) {
    msg <- paste('info: ', paste(list(...), collapse=''), '\n', sep='')
    cat(msg)
  }
}

warn <- function(...) {
  if (is.null(.log.info$level) || .log.info$level >= str2int('warn')) {
    msg <- paste('warning: ', paste(list(...), collapse=''), '\n', sep='')
    cat(msg)
  }
}

error <- function(...) {
  if (is.null(.log.info$level) || .log.info$level >= str2int('error')) {
    msg <- paste('error: ', paste(list(...), collapse=''), '\n', sep='')
    cat(msg)
  }
}
