grule <- function(...) {
  rule = eval(substitute(alist(...)))
  class(rule) = "GERule"
  return (rule)
}

gsrule <- function(...) {
  rule = list(...)
  class(rule) = c("GERule", "GEStringRule")
  return (rule)
}

gvrule <- function(vec) {
# create GERule for a vector
  do.call(grule, as.list(vec))
}

print.GERule <- function(x, ...) {
  i = 0
  for (item in x) {
    str = do.call(paste, as.list(as.character(GERule.EscapeDot(item))))
    cat(paste0('Rule ', i, ': '))
    cat(str)
    cat('\n')
    i = i + 1
  }
}

GERule.EscapeDot <- function(rule) {
  # convert
  # escape dot function for a given expression
  # .(x,x,x) -> "x,x,x"
  
  if (length(rule) > 1) {
    if (as.character(rule[[1]]) == '.') {
      str = as.character(as.expression(rule))
      return (substr(str, 3, nchar(str) - 1))
    } 
  }
  
  deparse(rule, width.cutoff=400)
}
