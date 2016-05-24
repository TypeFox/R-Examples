"[.mcmc" <- function (x, i, j, drop = missing(i)) 
{
  ## In S-PLUS the code is altered so that the user can
  ## pick out particular parameters by calling mcmc.obj[,c("param1", "param2")]
  xstart <- start(x)
  xthin <- thin(x)
  if (is.R()) {
    y <- NextMethod("[")
  }
  else {
    y <- as.matrix(x)[i,j]
  }
  if (length(y) == 0 || is.null(y)) 
    return(y)
  if (missing(i)) 
    return(mcmc(y, start = xstart, thin = xthin))
  else
    return(y)
}

"as.mcmc" <- function (x, ...) 
  UseMethod("as.mcmc")

"as.mcmc.default" <- function (x, ...) 
  if (is.mcmc(x)) x else mcmc(x)

"as.ts.mcmc" <- function (x, ...) 
{
  x <- as.mcmc(x)
  y <- ts(x, start = start(x), end = end(x), deltat = thin(x))
  attr(y, "mcpar") <- NULL
  return(y)
}

"start.mcmc" <- function (x, ...) 
{
  mcpar(as.mcmc(x))[1]
}

"end.mcmc" <- function (x, ...) 
{
  mcpar(as.mcmc(x))[2]
}

"frequency.mcmc" <- function (x, ...) 
{
  1/thin.mcmc(x)
}

"thin.mcmc" <- function (x, ...) 
{
  mcpar(as.mcmc(x))[3]
}

"is.mcmc" <- function (x) 
{
  if (inherits(x, "mcmc")) 
    if (length(dim(x)) == 3) 
      stop("Obsolete mcmc object\nUpdate with a command like\nx <- mcmcUpgrade(x)")
    else TRUE
  else FALSE
}

"mcmc" <- function (data = NA, start = 1, end = numeric(0), thin = 1) 
{
  if (is.matrix(data)) {
    niter <- nrow(data)
    nvar <- ncol(data)
  }
  else if (is.data.frame(data)) {
      if (!all(sapply(data, is.numeric))) {
         stop ("Data frame contains non-numeric values")
      }
      data <- as.matrix(data)
      niter <- nrow(data)
      nvar <- ncol(data)
  }
  else {
    niter <- length(data)
    nvar <- 1
  }
  thin <- round(thin)
  if (length(start) > 1) 
    stop("Invalid start")
  if (length(end) > 1) 
    stop("Invalid end")
  if (length(thin) != 1) 
    stop("Invalid thin")
  if (missing(end)) 
    end <- start + (niter - 1) * thin
  else if (missing(start)) 
    start <- end - (niter - 1) * thin
  nobs <- floor((end - start)/thin + 1.0) ### patch
  if (niter < nobs) 
    stop("Start, end and thin incompatible with data")
  else {
    end <- start + thin * (nobs - 1)
    if (nobs < niter) 
      data <- data[1:nobs, , drop = FALSE]
  }
  attr(data, "mcpar") <- c(start, end, thin)
  attr(data, "class") <- "mcmc"
  data
}

"print.mcmc" <- function (x, ...) 
{
  x.orig <- x
  cat("Markov Chain Monte Carlo (MCMC) output:\nStart =", start(x), 
      "\nEnd =", end(x), "\nThinning interval =", thin(x), "\n")
  attr(x, "mcpar") <- NULL
  attr(x, "class") <- NULL
  NextMethod("print", ...)
  invisible(x.orig)
}


"as.matrix.mcmc" <- function (x, iters = FALSE, ...) 
{
  y <- matrix(nrow = niter(x), ncol = nvar(x) + iters)
  var.cols <- iters + 1:nvar(x)
  if (iters) 
    y[, 1] <- as.vector(time(x))
  y[, var.cols] <- x
  rownames <- character(ncol(y))
  if (iters) 
    rownames[1] <- "ITER"
  rownames[var.cols] <- varnames(x, allow.null = FALSE)
  dimnames(y) <- list(NULL, rownames)
  return(y)
}

"time.mcmc" <- function (x, ...) 
{
  x <- as.mcmc(x)
  ts(seq(from = start(x), to = end(x), by = thin(x)), start = start(x), 
     end = end(x), deltat = thin(x))
}

"window.mcmc" <- function (x, start, end, thin, ...) 
{
  ts.eps <- getOption("ts.eps")
  xmcpar <- mcpar(x)
  xstart <- xmcpar[1]
  xend <- xmcpar[2]
  xthin <- xmcpar[3]
  if (missing(thin)) 
    thin <- xthin
  else if (thin%%xthin != 0) {
    thin <- xthin
    warning("Thin value not changed")
  }
  xtime <- as.vector(time(x))
  if (missing(start)) 
    start <- xstart
  else if (length(start) != 1) 
    stop("bad value for start")
  else if (start < xstart) {
    start <- xstart
    warning("start value not changed")
  }
  if (missing(end)) 
    end <- xend
  else if (length(end) != 1) 
    stop("bad value for end")
  else if (end > xend) {
    end <- xend
    warning("end value not changed")
  }
  if (start > end) 
    stop("start cannot be after end")
  if (all(abs(xtime - start) > abs(start) * ts.eps)) {
    start <- xtime[(xtime > start) & ((start + xthin) > xtime)]
  }
  if (all(abs(end - xtime) > abs(end) * ts.eps)) {
    end <- xtime[(xtime < end) & ((end - xthin) < xtime)]
  }
  use <- 1:niter(x)
  use <- use[use >= trunc((start - xstart)/xthin + 1.5) &
             use <= trunc((end - xstart)/xthin + 1.5) &
             (use - trunc((start- xstart)/xthin + 1.5))%%(thin%/%xthin) == 0]
  y <- if (is.matrix(x)) 
    x[use, , drop = FALSE]
  else x[use]
  return(mcmc(y, start=start, end=end, thin=thin))
}

"mcpar" <- function (x) 
{
  attr(x, "mcpar")
}

"mcmcUpgrade" <- function (x) 
{
  if (inherits(x, "mcmc")) {
    if (length(dim(x)) == 3) {
      nchain <- dim(x)[3]
      xtspar <- attr(x, "tspar")
      xstart <- xtspar[1]
      xend <- xtspar[2]
      xthin <- xtspar[3]
      out <- vector("list", nchain)
      for (i in 1:nchain) {
        y <- unclass(x)[, , 1, drop = TRUE]
        attr(y, "title") <- NULL
        attr(y, "tspar") <- NULL
        out[[i]] <- mcmc(y, start = xstart, end = xend, 
                         thin = xthin)
      }
      if (nchain == 1) 
        return(out[[1]])
      else return(mcmc.list(out))
    }
    else return(x)
  }
  else stop("Can't upgrade")
}

"thin" <-
function (x, ...)
  UseMethod("thin")

"set.mfrow" <-
function (Nchains = 1, Nparms = 1, nplots = 1, sepplot = FALSE) 
{
  ## Set up dimensions of graphics window: 
  ## If only density plots OR trace plots are requested, dimensions are: 
  ##	1 x 1	if Nparms = 1 
  ##	1 X 2 	if Nparms = 2 
  ##	2 X 2 	if Nparms = 3 or 4 
  ##	3 X 2 	if Nparms = 5 or 6 or 10 - 12 
  ##	3 X 3 	if Nparms = 7 - 9 or >= 13 
  ## If both density plots AND trace plots are requested, dimensions are: 
  ##	1 x 2	if Nparms = 1 
  ##	2 X 2 	if Nparms = 2 
  ##	3 X 2 	if Nparms = 3, 5, 6, 10, 11, or 12 
  ##	4 x 2	if Nparms otherwise 
  ## If separate plots are requested for each chain, dimensions are: 
  ##	1 x 2	if Nparms = 1 & Nchains = 2 
  ##	2 X 2 	if Nparms = 2 & Nchains = 2 OR Nparms = 1 & Nchains = 3 or 4 
  ##	3 x 2	if Nparms = 3 or >= 5 & Nchains = 2  
  ##		   OR Nchains = 5 or 6 or 10 - 12 (and any Nparms) 
  ##	2 x 3	if Nparms = 2 or 4 & Nchains = 3 
  ##	4 x 2   if Nparms = 4 & Nchains = 2  
  ##		   OR Nchains = 4 & Nparms > 1 
  ##	3 x 3	if Nparms = 3 or >= 5  & Nchains = 3  
  ##		   OR Nchains = 7 - 9 or >= 13 (and any Nparms)
  mfrow <- if (sepplot && Nchains > 1 && nplots == 1) {
    ## Separate plots per chain
    ## Only one plot per variable
    if (Nchains == 2) {
      switch(min(Nparms, 5),
             c(1,2),
             c(2,2),
             c(3,2),
             c(4,2),
             c(3,2))
    }
    else if (Nchains == 3) {
      switch(min(Nparms, 5),
             c(2,2),
             c(2,3),
             c(3,3),
             c(2,3),
             c(3,3))
    }
    else if (Nchains == 4) {
      if (Nparms == 1)
        c(2,2)
      else
        c(4,2)
    }
    else if (any(Nchains == c(5,6,10,11,12)))
      c(3,2)
    else if (any(Nchains == c(7,8,9)) || Nchains >=13)
      c(3,3)
      
  }
  else {
    if (nplots==1) {
      ## One plot per variable
      mfrow <- switch(min(Nparms,13),
                      c(1,1),
                      c(1,2),
                      c(2,2),
                      c(2,2),
                      c(3,2),
                      c(3,2),
                      c(3,3),
                      c(3,3),
                      c(3,3),
                      c(3,2),
                      c(3,2),
                      c(3,2),
                      c(3,3))
    }
    else {
      ## Two plot per variable
      ##
      mfrow <- switch(min(Nparms, 13),
                      c(1,2),
                      c(2,2),
                      c(3,2),
                      c(4,2),
                      c(3,2),
                      c(3,2),
                      c(4,2),
                      c(4,2),
                      c(4,2),
                      c(3,2),
                      c(3,2),
                      c(3,2),
                      c(4,2))
    }
  }
  return(mfrow)
}

head.mcmc <- function(x, n = 6L, ...) {
    window.mcmc(x, end=min(start.mcmc(x) + n * thin.mcmc(x), end.mcmc(x)))
}

tail.mcmc <- function(x, n = 6L, ...) {
    window.mcmc(x, start=max(end.mcmc(x) - n * thin.mcmc(x), start.mcmc(x)))
}
