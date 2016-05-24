get.unique <-
function(x, digits = 4L)
{
  check <- is.matrix(x)
  if(check)
    n <- nrow(x)
  else
    n <- length(x)
  if(mode(digits) == "numeric")
    xr <- round(x, digits = digits)
  else
    xr <- x
  iind <- !duplicated(xr)
  if(check) {
    xu <- xr[iind,]
    lu <- nrow(xu)
  } else {
    xu <- xr[iind]
    lu <- length(xu)
  }
  get <- getuit(xr, xu, n, lu)
  index <- get$ind

  return(list(x = xr, xu = xu, ind = index, iind = iind))
}

