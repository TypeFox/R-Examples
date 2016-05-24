
rvmatch <- function (x, table, nomatch=NA_integer_, incomparables=NULL) {
  rvmapply(base::match, x=x, table=table, nomatch=nomatch, MoreArgs=list(incomparables=incomparables))
}

"%*in*%" <- function (x, y) {
  if (! is.rv(x) && ! is.rv(y)) {
    return(.Primitive("%in%")(x, y))
  }
  z <- rvmatch(x, table=y, nomatch=0L)
  return(z > 0)
}
