"interpp.old"<-function(x, y, z, xo, yo, ncp = 0, extrap = FALSE,
                    duplicate = "error", dupfun = NULL)
{
  if(!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
    stop("missing values and Infs not allowed")
  if(is.null(xo))
    stop("xo missing")
  if(is.null(yo))
    stop("yo missing")
  if(ncp>25){
    ncp <- 25
    cat("ncp too large, using ncp=25\n")
  }
  drx <- diff(range(x))
  dry <- diff(range(y))
  if(drx == 0 || dry == 0)
    stop("all data collinear")	# other cases caught in Fortran code
  if(drx/dry > 10000 || drx/dry < 0.0001)
    stop("scales of x and y are too dissimilar")
  n <- length(x)
  np <- length(xo)
  if(length(yo)!=np)
    stop("length of xo and yo differ")
  if(length(y) != n || length(z) != n)
    stop("Lengths of x, y, and z do not match")
  xy <- paste(x, y, sep =",")
  i <- match(xy, xy)
  if(duplicate=="user" && !is.function(dupfun))
    stop("duplicate=\"user\" requires dupfun to be set to a function")
  if(duplicate!="error")
    {
      centre <- function(x) {
        switch(duplicate,
               mean = mean(x),
               median = median(x),
               user = dupfun(x)
               )
      }
      if(duplicate!="strip"){
        z <- unlist(lapply(split(z,i), centre))
        ord <- !duplicated(xy)
        x <- x[ord]
        y <- y[ord]
        n <- length(x)
      }
      else{
        ord <- (hist(i,plot=FALSE,freq=TRUE,breaks=seq(0.5,max(i)+0.5,1))$counts==1)
        x <- x[ord]
        y <- y[ord]
        z <- z[ord]
        n <- length(x)
      }
    }
  else
    if(any(duplicated(xy)))
      stop("duplicate data points")
  zo <- rep(0, np)
  storage.mode(zo) <- "double"
  miss <- !extrap	#if not extrapolating use missing values
  misso <- seq(miss, np)
  if(extrap & ncp == 0)
    warning("Cannot extrapolate with linear option")
  ans <- .Fortran("idbvip",
                  as.integer(1),
                  as.integer(ncp),
                  as.integer(n),
                  as.double(x),
                  as.double(y),
                  as.double(z),
                  as.integer(np),
                  x = as.double(xo),
                  y = as.double(yo),
                  z = zo,
                  integer((31 + ncp) * n + np),
                  double(8 * n),
                  misso = as.logical(misso),
                  PACKAGE = "akima")
  temp <- ans[c("x", "y", "z", "misso")]
  temp$z[temp$misso]<-NA
  temp[c("x", "y", "z")]
}
