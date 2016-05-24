# Parse and evaluate an RSurvey expression.

EvalFunction <- function(txt, cols) {
  d <- list()
  ids <- vapply(cols, function(i) i$id, "")
  for (i in seq_along(ids)) {
    if (regexpr(paste0("\"", ids[i], "\""), txt, fixed=TRUE)[1] >= 0) {
      if (is.na(cols[[i]]$index))
        d[[i]] <- EvalFunction(cols[[i]]$fun, cols)
      else
        d[[i]] <- Data("data.raw")[[cols[[i]]$index]]
    }
  }
  fun <- txt
  for (i in seq_along(ids))
    fun <- gsub(paste0("\"", ids[i], "\""), paste0("DATA[[", i, "]]"), fun,
                fixed=TRUE)
  fun <- eval(parse(text=paste0("function(DATA) {", fun, "}")))
  ans <- try(fun(d), silent=TRUE)
  return(ans)
}
