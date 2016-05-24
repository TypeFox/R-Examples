
lambda <- function(t, x, y, param, object)
{
  if (!is.numeric(t) || !is.numeric(x) || !is.numeric(y))
    stop(paste("Arguments", sQuote(t), ",", sQuote(x), "and", sQuote(y), "must be numeric."))
  if (length(t) != length(x) || length(t) != length(y) || length(x) != length(y))
    stop(paste("Arguments", sQuote(t), ",", sQuote(x), "and", sQuote(y), "must be of the same length."))

  spatstat::verifyclass(object, "catalog")
  revents <-  object$revents

  storage.mode(revents) <- "double"
  theta <- sqrt(param)
  out <- numeric(length(t))
  for (i in 1:length(t))
  {
    if (t[i] < revents[1,1])
      out[i] <- 0
    else
    {
      out[i] <- .Call("lambdax", as.double(t[i]), as.double(x[i]), as.double(y[i]),
                 as.double(theta), revents, PACKAGE="ETAS")[[1]]
    }
  }
  return(out)
}

