
"[.rv" <- function (x, ..., drop = TRUE)
{
  cx <- class(x)
  X <- NextMethod()
  class(X) <- cx
  return(X)
}

"[.rvsummary" <- function (x, ..., drop = TRUE)
{
  q <- attr(x, "quantiles")
  cx <- class(x)
  x <- NextMethod()
  class(x) <- cx
  attr(x, "quantiles") <- q
  return(x)
}


"[<-.rvsummary" <- function (x, ..., value = NULL)
{
  cx <- class(x)
  q <- attr(x, "quantiles")
  value <- as.rvsummary(value, quantiles=q)
  X <- .Primitive("[<-")(unclass(x), ..., value=value)
  class(X) <- cx
  return(X)
}

"[<-.rv" <- function (x, ..., value = NULL)
{
  cx <- class(x)
  value <- as.rvobj(value)
  X <- .Primitive("[<-")(unclass(x), ..., value=value)
  class(X) <- cx
  return(X)
}

"impute<-" <- function(x, ..., value) {
  if (! is.rvobj(x) && ! is.rvobj(value)) {
    x[...] <- value
  } else if (is.rvsummary(x) || is.rvsummary(value)) {
    x <- as.rvsummary(x)
    value <- as.rvsummary(value)
    x[...] <- value
  } else if (is.rv(x) || is.rv(value)) {
    x <- as.rv(x)
    value <- as.rv(value)
    x[...] <- value
  }
  return(x)
}





