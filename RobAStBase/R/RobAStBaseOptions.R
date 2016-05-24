.RobAStBaseOptions <- list(
    kStepUseLast = FALSE,
    withUpdateInKer = FALSE,
    IC.UpdateInKer = NULL,
    all.verbose = FALSE,
    withICList = FALSE,
    withPICList = FALSE
)

RobAStBaseOptions <- function(...) {
  if (nargs() == 0) return(.RobAStBaseOptions)
  current <- .RobAStBaseOptions
  temp <- list(...)
  if (length(temp) == 1 && is.null(names(temp))) {
    arg <- temp[[1]]
    switch(mode(arg),
           list = temp <- arg,
           character = return(.RobAStBaseOptions[arg]),
           stop("invalid argument: ", sQuote(arg)))
  }
  if (length(temp) == 0) return(current)
  n <- names(temp)
  if (is.null(n)) stop("options must be given by name")
  changed <- current[n]
  current[n] <- temp
  env <- if (sys.parent() == 0) asNamespace("RobAStBase") else parent.frame()
  assign(".RobAStBaseOptions", current, envir = env)
  invisible(current)
}

getRobAStBaseOption <- function(x) RobAStBaseOptions(x)[[1]]
