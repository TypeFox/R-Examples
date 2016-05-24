match.bitarg <- function(arg, choices) {
  if (missing(choices)) {
    formal.args <- formals(sys.function(sys.parent()))
    choices <- eval(formal.args[[deparse(substitute(arg))]])
  }
  resp <- logical(length(choices))
  names(resp) <- choices

  if (eval.parent(substitute(missing(arg))) || is.null(arg)) {
    resp[1L] <- TRUE
    return(resp)
  } else if (!is.character(arg))
    stop("'arg' must be NULL or a character vector")
  if (length(arg) == 0L)
    stop("'arg' must be of length  >= 1")
  i <- unique(pmatch(arg, choices, nomatch=0L, duplicates.ok=TRUE))
  if (all(i == 0L))
    stop(gettextf("'arg' should be one of %s", paste(dQuote(choices),
                                                     collapse = ", ")),
         domain = NA)
  
  i <- i[i > 0L]
  resp[i] <- TRUE

  return(resp)
}
