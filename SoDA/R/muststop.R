muststop <- function(expr, silent = TRUE) {
    tryExpr <- substitute(tryCatch(expr, error=function(cond)cond))
    value <- eval.parent(tryExpr)
    if(inherits(value, "error")) {
        if(!silent)
          message("muststop reports: ", value)
        invisible(value)
    }
    else
      stop(gettextf("The expression  %s should have thrown an error, but instead returned an object of class \"%s\"",
           deparse(substitute(expr))[[1]], class(value)))
}
