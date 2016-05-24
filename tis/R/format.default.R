## functions in the tis namespace will find this 
format.default <- function(x, ...){
  ## BUGFIX: don't lose attributes of x
  res <- get("format.default", pos = "package:base")(x, ...)
  attributes(res) <- attributes(x)
  res
}
