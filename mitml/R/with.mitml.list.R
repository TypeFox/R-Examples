with.mitml.list <- function(data, expr, ...){
# evaluates an expression for a list of data sets

  expr <- substitute(expr)
  parent <- parent.frame()

  out <- lapply(data, function(x) eval(expr, x, parent))
  class(out) <- c("mitml.result","list")
  out

}
