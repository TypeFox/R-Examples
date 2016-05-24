.distrModOptions <- list(
                         show.details = "maximal"
                         )


distrModOptions <- function(...) {
  if (nargs() == 0) return(.distrModOptions)
  current <- .distrModOptions
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.distrModOptions[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- temp
  env <- if (sys.parent() == 0) asNamespace("distrMod") else parent.frame()
  assign(".distrModOptions", current, envir = env)
  invisible(current)
}

getdistrModOption <- function(x) distrModOptions(x)[[1]]
distrModoptions <- distrModOptions
