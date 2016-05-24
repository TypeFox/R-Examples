
unlist.rv <- function (x, recursive = TRUE, use.names = TRUE) {
  ## NAME
  ##   unlist.rv - Flatten Lists That (May) Contain Random Variables 
  y <- NULL
  ix <- seq(along=x)
  xn <- names(x)
  .paste <- function (name, x) {
     nbrs <- .dim.index(x, leftadjust=FALSE)
     paste(name, nbrs, names(x), sep="")
  }
  if (recursive) {
    for (i in ix) {
      nx <- xn[i]
      if (use.names && is.null(nx)) nx <- "."
      if (!is.rv(x[[i]]) && is.list(x[[i]])) {
        new.y <- unlist.rv(x[[i]], recursive=TRUE, use.names=use.names)
      } else {
        new.y <- x[[i]]
      }
      if (is.null(names(new.y))) {
        new.names <- .paste(nx, new.y)
      } else {
        new.names <- paste(nx, ".", names(new.y), sep="")
      }
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  } else {
    for (i in ix) {
      nx <- xn[i]
      new.y <- x[[i]]
      new.names <- .paste(nx, new.y)
      yn <- names(y)
      y <- c(y, new.y)
      names(y) <- c(yn, new.names)
    }
  }
  return(y)
}
