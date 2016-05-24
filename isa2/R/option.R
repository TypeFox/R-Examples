
if (!exists("isa.options")) { isa.options <- new.env() }

isa.option <- function(...) {
  args <- list(...)
  nam <- names(args)

  if (length(args) == 0) {
    return(as.list(isa.options))
  }
  
  if (is.null(nam)) {
    ## query
    if (length(args) != 1) {
      warning("Can't query many options at once, ignoring the rest")
    }
    if (!is.character(args[[1]]) || length(args[[1]]) != 1) {
      stop("Expected character of length one")
    }
    return(isa.options[[ args[[1]] ]])
  } else {
    ## set, everything must be named
    if (any(nam=="")) {
      stop("Some options are not named")
    }
    for (i in seq_along(nam)) {
      isa.options[[ nam[i] ]] <- args[[i]]
    }
    return(invisible(isa.options))
  }
}
