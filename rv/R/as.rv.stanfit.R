
as.rv.stanfit <- function (x, list.=TRUE, ...) {
  f <- get("extract", mode="function")
  L <- lapply(f(x), rvsims)
  P <- x@par_dims
  for (name in names(L)) {
    x <- L[[name]]
    n <- P[[name]]
    if (length(n) == 0) {
      names(x) <- name
    } else if (length(n) == 1) {
      names(x) <- paste(name, "[", 1:n, "]", sep="")
    } else {
      C <- lapply(n, function (k) 1:k)
      A <- as.list(do.call(expand.grid, C))
      A$sep <- ","
      indices <- do.call(paste, A)
      names(x) <- paste(name, "[", indices, "]", sep="")
    }
    L[[name]] <- x
  }
  if (! list.) {
    L <- do.call(c.rv, L)
  }
  return(L)
}

