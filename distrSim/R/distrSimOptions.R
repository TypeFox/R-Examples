distrSimoptions <- function(...) {
  if (nargs() == 0) return(.distrSimoptions)
  current <- .distrSimoptions
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.distrSimoptions[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- temp
  if (sys.parent() == 0) env <- asNamespace("distrSim") else env <- parent.frame()
  assign(".distrSimoptions", current, envir = env)
  invisible(current)
}

getdistrSimOption<-function(x)
distrSimoptions(x)[[1]]
