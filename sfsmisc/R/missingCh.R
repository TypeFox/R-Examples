missingCh <- function(x, envir = parent.frame()) {
    stopifnot(is.character(x))
    eval(substitute(missing(VAR), list(VAR=as.name(x))),
         envir = envir)
}
