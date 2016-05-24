swap <- function(x, y){
  x.sub <- substitute(x)
  y.sub <- substitute(y)
  x.val <- x
  e <- parent.frame()
  do.call("<-", list(x.sub, y.sub), envir = e)
  do.call("<-", list(y.sub, x.val), envir = e)
}

