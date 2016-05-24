fmatch <- function(x, table, nomatch = NA_integer_, incomparables = NULL)
  .Call("fmatch", x, table, nomatch, incomparables, PACKAGE = "fastmatch")
