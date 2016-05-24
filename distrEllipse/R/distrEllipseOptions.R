distrEllipseoptions <- function(...) {
  if (nargs() == 0) return(.distrEllipseoptions)
  current <- .distrEllipseoptions
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.distrEllipseoptions[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- temp
  if (sys.parent() == 0) env <- asNamespace("distrEllipse") else env <- parent.frame()
  assign(".distrEllipseoptions", current, envir = env)
  invisible(current)
}

getdistrEllipseOption<-function(x)
distrEllipseoptions(x)[[1]]
