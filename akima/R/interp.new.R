interp.new <-
  function(x, y, z,
           xo = seq(min(x), max(x), length = 40),
           yo = seq(min(y), max(y), length = 40), linear = FALSE,
           ncp = NULL, extrap = FALSE, duplicate = "error", dupfun = NULL)
{
  if(!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
    stop("missing values and Infs not allowed")
  if(!is.null(ncp)) {
    if(ncp != 0)
	warning("'ncp' not supported, it is automatically choosen by Fortran code\n")
    else
	linear <- TRUE
  }
  if(linear)
    stop("linear interpolation not implemented in interp.new().\n",
	 "use 'interp()' (or 'interp.old()').")

  drx <- diff(range(x))
  dry <- diff(range(y))
  if(drx == 0 || dry == 0)
    stop("all data collinear")    # other cases caught in Fortran code
  if(drx/dry > 10000 || drx/dry < 0.0001)
    stop("scales of x and y are too dissimilar")
  n <- length(x)
  nx <- length(xo)
  ny <- length(yo)
  if(length(y) != n || length(z) != n)
    stop("Lengths of x, y, and z do not match")

  xy <- paste(x, y, sep = ",")# trick for 'duplicated' (x,y)-pairs
  if(duplicate == "error") {
    if(any(duplicated(xy)))
      stop("duplicate data points: need to set 'duplicate = ..' ")
  }
  else { ## duplicate != "error"

    i <- match(xy, xy)
    if(duplicate == "user")
      dupfun <- match.fun(dupfun)#> error if it fails

    ord <- !duplicated(xy)
    if(duplicate != "strip") {
      centre <- function(x)
        switch(duplicate,
               mean = mean(x),
               median = median(x),
               user = dupfun(x))
      z <- unlist(lapply(split(z,i), centre))
    } else {
      z <- z[ord]
    }
    x <- x[ord]
    y <- y[ord]
    n <- length(x)
  }

  zo <- matrix(0, nx, ny)
  storage.mode(zo) <- "double"
  miss <- !extrap             # if not extrapolating, set missing values
  misso <- matrix(TRUE, nx, ny)# hmm, or rather 'miss' ??

  if(extrap && if(is.null(ncp)) linear else (ncp == 0))
      warning("Cannot extrapolate with linear option")

  ans <- .Fortran("sdsf3p",
                  as.integer(1),
                  as.integer(n),
                  as.double(x),
                  as.double(y),
                  as.double(z),
                  as.integer(nx),
                  x = as.double(xo),
                  as.integer(ny),
                  y = as.double(yo),
                  z = zo,
                  ier = integer(1),
                  double(36 * n),
                  integer(25 * n),
                  extrap = as.logical(misso),
                  near = integer(n),
                  nxt = integer(n),
                  dist = double(n),
                  PACKAGE = "akima")[c("x", "y", "z", "extrap")]
  if(miss)
    ans$z[ans$extrap] <- NA

  ans[c("x", "y", "z")]
}
