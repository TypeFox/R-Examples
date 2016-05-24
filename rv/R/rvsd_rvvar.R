

rvvar <- function (x)
{
  UseMethod("rvvar")
}

rvvar.rv <- function (x) # NOEXPORT
{
  S <- sims(x)
  m <- colMeans(S, na.rm=TRUE)
  ns <- rvnsims(x)
  v <- ((colSums(S^2)-ns*(m^2))/(ns-1))
  v[ns==1] <- 0
  names(v) <- names(x)
  dim(v) <- dim(x)
  dimnames(v) <- dimnames(x)
  return(v)
}

rvvar.rvsummary <- function (x) # NOEXPORT
{
  return(unlist(rvattr(x, "sd"), use.names=TRUE)^2)
}

rvvar.default <- function (x) # NOEXPORT
{
  rep.int(0, length(x))
}

rvsd <- function (x)
{
  UseMethod("rvsd")
}

rvsd.rv <- function (x) # NOEXPORT
{
  sqrt(rvvar(x))
}

rvsd.rvsummary <- function (x) # NOEXPORT
{
  unlist(rvattr(x, "sd"), use.names=TRUE)
}

rvsd.default <- function (x) # NOEXPORT
{
  rep.int(0, length(x))
}



