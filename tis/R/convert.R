convert <- function(x, tif, method = "constant", observed. = observed(x), 
                    basis. = basis(x), ignore = F){
  ## converts series x to the ti frequency tif by one of the
  ## methods outlined in FAME documentation
  wasMatrix <- is.matrix(x)
  xdm <- dimnames(x)
  x <- as.matrix(as.tis(x))
  n <- nrow(x)
  m <- ncol(x)
  idates <- ti(x)

  method. <- match.arg(method, c("discrete", "constant", "linear", "cubic"))

  ## make sure basis. and observed. are set to proper values
  if(is.null(observed.)) observed. <- "averaged"
  observed(x) <- observed.
  observed.   <- observed(x)
  if(observed. %in% c("high", "low"))
    method. <- "discrete"
  
  if(is.null(basis.))    basis.    <- "daily"
  basis(x) <- basis.
  basis.   <- basis(x)

  bTif <- basis(x)
  xTif <- tif(x)
  agg <- (frequency(x) >= frequency(ti(tif = tif)))
  
  if((tif == xTif) || (!agg && xTif == bTif)) return(x)

  convert.discrete <- function(x, tif, bTif, ignore = F){
    if(agg){
      odates <- ti(jul(idates), tif)
      agg.fun <- switch(observed.,
                        beginning  = function(v, arg2) v[1],
                        ending     = function(v, arg2) v[length(v)],
                        high       = function(v, arg2) max( v, na.rm = arg2),
                        low        = function(v, arg2) min( v, na.rm = arg2),
                        summed     = function(v, arg2) sum( v, na.rm = arg2),
                        averaged   = ,
                        annualized = function(v, arg2) mean(v, na.rm = arg2))
      
      y <- tis(tapply(unclass(x),
                      list(rep(unclass(odates), ncol(x)), col(x)),
                      agg.fun,
                      arg2 = ignore),
               start = odates[1]) 
      
      if(!ignore){ ## chop months with incomplete data
        ostart <- asTi(ifelse(observed. == "ending",
                               start(y),
                               ti(jul(start(x) - 1), tif) + 1))
        oend <- asTi(ifelse(observed. == "beginning", end(y),
                             oend <- ti(jul(end(x) + 1), tif) - 1))
      
        if(ostart > oend) stop("Not enough observations to fill in an output period")
        y <- window(y, start = ostart, end = oend)
      }
    }
    else {
      ## for simplicity, pretend that input is monthly and output weekly
      ## First and Lasts Weeks that End in the Months
      fwem <- ti(jul(idates - 1) + 1, tif)
      lwem <- ti(jul(idates) + 1, tif) - 1
    
      switch(observed.,
             beginning = {
               y <- tis(matrix(NA, ncol = ncol(x), nrow = 1 + max(fwem) - min(fwem)),
                        start = min(fwem))
               y[fwem,] <- x
             },
             ending = {
               y <- tis(matrix(NA, ncol = ncol(x), nrow = 1 + max(lwem) - min(lwem)),
                        start = min(lwem))
               y[lwem,] <- x
             },
             summed = {
               runs <- rep(lwem - fwem + 1, ncol(x))
               y <- tis(matrix(NA, ncol = ncol(x), nrow = 1 + max(lwem) - min(fwem)),
                        start = min(fwem))
               y[1:nrow(y),] <- rep(unclass(x)/runs, runs)
             },
             averaged = ,
             high = ,
             low = ,
             annualized = {
               runs <- rep(lwem - fwem + 1, ncol(x))
               y <- tis(matrix(NA, ncol = ncol(x), nrow = 1 + max(lwem) - min(fwem)),
                        start = min(fwem))
               y[1:nrow(y),] <- rep(unclass(x), runs)
             })
    }
    y
  }

  convert.constant <- function(x, tif, bTif, ignore = F){
    switch(observed.,
           beginning = {
             ibasis <- ti(jul(idates - 1), bTif) + 1
           
             if(ignore) ## only need 1 datum for 1st output period
               ostart <- ti(jul(ibasis[1]), tif)
             else ## need complete data for 1st output period
               ostart <- ti(jul(ibasis[1] - 1), tif) + 1
           
             oend   <- ti(jul(ibasis[n]), tif)
             odates <- seq(ostart, oend)
             obasis <- ti(jul(odates - 1), bTif) + 1
             vals <- matrix(0, length(obasis), m)
             for(i in 1:m){
               vals[,i] <- approx(x = ibasis, 
                                  y = x[,i], 
                                  xout = obasis,
                                  method = "constant",
                                  f = 0,
                                  rule = 2)$y
             }
             yy <- tis(vals, start=ostart)
           },
           ending = {
             ibasis <- ti(jul(idates), bTif)
           
             ostart <- ti(jul(idates[1] - 1) + 1, tif)
           
             if(ignore) ## only need 1 datum for last output period
               oend <- ti(jul(idates[n]), tif)
             else ## need complete data for last output period
               oend <- ti(jul(idates[n]) + 1, tif) - 1
           
             odates <- seq(ostart, oend)
             obasis <- ti(jul(odates), bTif)
             vals <- matrix(0, length(obasis), m)
             for(i in 1:m){
               vals[,i] <- approx(x = ibasis, 
                                  y = x[,i], 
                                  xout = obasis,
                                  method = "constant",
                                  f = 1,
                                  rule = 2)$y
             }
             yy <- tis(vals, start=ostart)
           },
           summed = , averaged =,
           annualized = {
             if(ignore){ ## only need one input datum for an output period
               ostart <- ti(jul(idates[1] - 1) + 1, tif)
               oend <- ti(jul(idates[n]), tif)
             }
             else { ## need all input datums for an output period
               ostart <- ti(jul(idates[1] - 1), tif) + 1
               oend <- ti(jul(idates[n]) + 1, tif) - 1
             }
             if(ostart > oend)
               stop("Not enough observations to fill in an output period")

             odates <- seq(ostart, oend)
           
             ## Basis-day borders for the input and output series
             iborders <- ti(jul(c(idates[1]-1, idates)), bTif)
             oborders <- ti(jul(c(odates[1]-1, odates)), bTif)
           
             idays <- diff(iborders) ## basis days per period of input series
             odays <- diff(oborders) ## basis days per period of output series
           
             if(observed. == "summed") x <- x/idays
           
             vals <- matrix(0, length(odays), m)
             for(i in 1:m){
               vals[,i] <- lintegrate(x = iborders,
                                      y = c(x[,i], 0),
                                      xint = oborders,
                                      stepfun = T,
                                      rule = 1)
             }
             if(observed. == "averaged" || observed. == "annualized") 
               vals <- vals/rep(odays, m)
             yy <- tis(vals, start = ostart)
           })
    yy
  }

  convert.linear <- function(x, tif, bTif, ignore = F){
    ## This is very similar to convert.constant
    switch(observed.,
           beginning = {
             ibasis <- ti(jul(idates - 1), bTif) + 1
           
             ostart <- ti(jul(ibasis[1] - 1), tif) + 1
             oend   <- ti(jul(ibasis[n]), tif)
           
             odates <- seq(ostart, oend)
             obasis <- ti(jul(odates - 1), bTif) + 1
             vals <- matrix(0, length(obasis), m)
             for(i in 1:m){
               vals[,i] <- approx(x = ibasis, 
                                  y = x[,i], 
                                  xout = obasis,
                                  rule = 2)$y
             }
             yy <- tis(vals, start=ostart)
           },
           ending = {
             ibasis <- ti(jul(idates), bTif)
           
             ostart <- ti(jul(idates[1]), tif)
             oend <- ti(jul(idates[n]) + 1, tif) - 1
           
             odates <- seq(ostart, oend)
             obasis <- ti(jul(odates), bTif)
             vals <- matrix(0, length(obasis), m)
             for(i in 1:m){
               vals[,i] <- approx(x = ibasis, 
                                  y = x[,i], 
                                  xout = obasis,
                                  rule = 2)$y}
             yy <- tis(vals, start=ostart)
           },
           summed = , averaged = ,
           annualized = {
             if(ignore){ ## only need one input datum for an output period
               ostart <- ti(jul(idates[1] - 1) + 1, tif)
               oend <- ti(jul(idates[n]), tif)
             }
             else { ## need all input datums for an output period
               ostart <- ti(jul(idates[1] - 1), tif) + 1
               oend <- ti(jul(idates[n]) + 1, tif) - 1
             }
             odates <- seq(ostart, oend)
           
             ## Basis-day borders for the input and output series
             iborders <- ti(jul(c(idates[1]-1, idates)), bTif)
             oborders <- ti(jul(c(odates[1]-1, odates)), bTif)
           
             idays <- diff(iborders) ## basis days per period of input series
             odays <- diff(oborders) ## basis days per period of output series
           
             if(observed. == "averaged" || observed. == "annualized") 
               x <- x*rep(idays, m)
           
             vals <- matrix(0, length(odays), m)
             xy <- ilspline(iborders, x)
             for(i in 1:m){
               vals[,i] <- lintegrate(x = xy$x,
                                      y = xy$y[,i],
                                      xint = oborders,
                                      rule = 1)
             }
             if(observed. == "averaged" || observed. == "annualized") 
               vals <- vals/rep(odays, m)
             yy <- tis(vals, start = ostart)
           })
    yy
  }

  convert.cubic <- function(x, tif, bTif, ignore = F){
    switch(observed.,
           beginning = {
             ibasis <- ti(jul(idates - 1), bTif) + 1
           
             ostart <- ti(jul(ibasis[1] - 1), tif) + 1
             oend   <- ti(jul(ibasis[n]), tif)
           
             odates <- seq(ostart, oend)
             obasis <- ti(jul(odates - 1), bTif) + 1
           
             ostart.b <- obasis[1]
             oend.b   <- obasis[length(obasis)]
           
             vals <- matrix(0, length(obasis), m)
             for(i in 1:m){
               basis.series <- tis(spline(x = ibasis,
                                          y = x[,i],
                                          n = oend.b - ostart.b + 1,
                                          xmin = ostart.b,
                                          xmax = oend.b)$y,
                                   start=ostart.b)
               vals[,i] <- basis.series[obasis]
             }
             yy <- tis(vals, start=ostart)
           },
           ending = {
             ibasis <- ti(jul(idates), bTif)
           
             ostart <- ti(jul(idates[1]), tif)
             oend <- ti(jul(idates[n]) + 1, tif) - 1
           
             odates <- seq(ostart, oend)
             obasis <- ti(jul(odates), bTif)
           
             ostart.b <- obasis[1]
             oend.b   <- obasis[length(obasis)]
           
             vals <- matrix(0, length(obasis), m)
             for(i in 1:m){
               basis.series <- tis(spline(x = ibasis,
                                          y = x[,i],
                                          n = oend.b - ostart.b + 1,
                                          xmin = ostart.b,
                                          xmax = oend.b)$y,
                                   start=ostart.b)
               vals[,i] <- basis.series[obasis]
             }
             yy <- tis(vals, start=ostart)
           },
           summed = , averaged = ,
           annualized = {
             stop(paste("cubic method not implemented for observed =",
                        observed.))
           })
    yy
  }
  y <- do.call(paste("convert", method., sep = "."),
               list(x = x, tif = tif, bTif = bTif, ignore = ignore))

  if(wasMatrix){
    if(is.null(xdm)) dimnames(y) <- NULL
    else             dimnames(y) <- list(character(0), xdm[[2]])
  }
  else {
    dimnames(y) <- NULL
    dim(y) <- NULL
  }
  
  observed(y) <- observed.
  basis(y) <- basis.
  class(y) <- "tis"
  return(y)
}

basis <- function(x) attr(x, "basis")

"basis<-" <- function(x, value){
  if(is.null(value))
    attr(x, "basis") <- NULL
  else attr(x, "basis") <- match.arg(value, c("daily", "business"))
  invisible(x)
}

observed <- function(x) attr(x, "observed")

"observed<-" <- function(x, value){
  if(is.null(value))
    attr(x, "observed") <- NULL
  else attr(x, "observed") <- match.arg(value,
                                        c("beginning", "ending", "high",
                                          "low","summed", "averaged", "annualized"))
  invisible(x)
}


ilspline <- function(xint, w){
  ## Given xint[0] < xint[1] < xint[2] < .... < xint[N]  
  ## and a matrix W with M columns (w1, w2, ..., wM)
  ## of length N [i.e., wi = (wi[1], wi[2], ...., wi[N])]
  ## returns the N-vector x, and NxM matrix Y with columns
  ## (y1, y2, ..., yM) such that
  ## (i)   x[j] = (xint[j-1] + xint[j])/2
  ## (ii)  the linear spline S that passes through the x,yi pairs 
  ##       (and is extended to xint[0] and xint[n] by linear extrapolation)
  ##       integrates over [xint[j-1], xint[j]] to wi[j]
  w <- as.matrix(w)
  n <- nrow(w)
  if(length(xint) != n + 1) stop("length(xint) != 1 + nrow(w)")
  z <- xint
  zstart <- if(is.tis(z)) start(z) else NULL
  wstart <- if(is.tis(w)) start(w) else NULL
  if(n > 400){
    y <- w*0
    x <- as.vector(z[-1]*0)
    splits <- seq(200, n-11, by=200)
    starts <- c(1, splits-9)
    ends <- c(splits + 10, n)
    rstarts <- c(0, splits) + 1
    rends <- c(splits, max(ends))
    istarts <- 1 + rstarts - starts
    iends <- 1 + rends - starts
    for(i in 1:(length(splits)+1)){
      wdex <- starts[i]:ends[i]
      zdex <- c(wdex, ends[i]+1)
      rdex <- rstarts[i]:rends[i]
      idex <- istarts[i]:iends[i]
      xyi <- ilspline(z[zdex], w[wdex,])
      y[rdex,] <- xyi$y[idex,]
      x[rdex] <- xyi$x[idex]
    }
  }
  else {
    z <- tis(z, start=ti(19010101, "daily"))
    d <- diff(z)
    q0 <- d/4
    l <- q0/(lag(q0, -1) + q0)
    r <- q0/(q0 + lag(q0))
    M <- diag(c(d[1], q0*(2*lag(r, -1) + l + r + 2*lag(l)), d[n]))
    rm <- row(M)
    cm <- col(M)
    M[rm == cm + 1] <- q0*l
    M[rm == cm - 1] <- q0*r
    M[1,2] <- M[n,n-1] <- 0
    class(M) <- "tridiag"
    y <- solve(M, w)
    x <- as.vector((z + lag(z, -1))/2)
  }
  if(!is.null(wstart)){
    y <- tis(y, start=wstart)
  }
  if(!is.null(zstart)){
    x <- tis(x, start=zstart+1)
  }
  list(x = x, y = y)
}

lintegrate <- function(x, y, xint, stepfun = F, rule = 0){
  ## Integrates the linear spline F defined by (x,y) over
  ## the xint intervals
  ## rule == 0  --> F(z) = 0 for z outside the range of x
  ## rule == NA --> F(z) = NA for z outside the range of x
  ## rule == 1  --> F(z) extended for z outside the range of x
  ## 
  ## If stepfun == T, then F(z) is assumed to be a left-continuous step function
  ## and the last value of y is never accessed.
  ##
  ## (x[i], y[i]) pairs with either x[i] == NA or y[i] == NA are ignored
  if(any(na.spots <- (is.na(x) | is.na(y)))){
	x <- x[!na.spots]
	y <- y[!na.spots]
  }
  xord <- order(x)
  x <- x[xord]
  y <- y[xord]
  xint <- as.double(sort(xint))
  xlen <- length(x)
  ilen <- length(xint)
  crule <- as.integer( ifelse( is.na(rule), 0, rule))
  z <- .C("lintegrate", 
		  x = as.double(x - xint[1]), 
		  y = as.double(y), 
		  xlen = xlen, 
		  integrals = xint - xint[1], 
		  ilen = ilen, 
		  stepfun = as.integer(stepfun),
		  rule = crule,
          PACKAGE = "tis")
  integrals <- z$integrals[-ilen]
  if(is.na(rule)){
    na.spots <- c((1:ilen)[xint < x[1]], (1:ilen)[xint > x[xlen]] - 1)
    integrals[na.spots] <- NA
  }
  integrals
}

