

as.list.rv <- function (x, ...) { ## 
  L <- vector(mode="list", length=length(x))
  for (i in seq_along(x)) {
    L[[i]] <- x[i]
  }
  return(L)
}
