sharpen.wmppp <- function(...)
{
  X <- sharpen.ppp(...)
  class(X) <- c("wmppp", "ppp")
  return(X)
}

split.wmppp <- function(...)
{
  Xlist <- split.ppp(...)
  as.wmppplist <- function(X) {
    class(X) <- c("wmppp", "ppp")
    return(X)
  }
  Xlist <- lapply(Xlist, as.wmppplist)
  return(Xlist)
}

superimpose.wmppp <- function(...)
{
  X <- superimpose.ppp(...)
  class(X) <- c("wmppp", "ppp")
  return(X)
}

unique.wmppp <- function(...)
{
  X <- unique.ppp(...)
  class(X) <- c("wmppp", "ppp")
  return(X)
}
