"[.mcmc.list" <- function (x, i, j, drop = TRUE) 
{
    ## In S-PLUS the code is altered so that the user can
    ## pick out particular parameters by calling
    ## mcmc.obj[,c("param1", "param2")]
  
    ## Trying to squeeze too much functionality in here
    ## x[p:q] will subset the list
    ## x[p,], x[,q], x[p,q] will be recursively applied to
    ## the elements of the list, even if they are vectors

    if (nargs() < 3 + !missing(drop)) {
        ## Subset the list
        if (is.R()) {
            y <- NextMethod("[")
        }
        else {
            y <- as.matrix(x)[i,j]
        }
    }
    else {
        ## Subset the elements of the list
        y <- vector("list", length(x))
        names(y) <- names(x)
        for (k in 1:length(y)) {
            y[[k]] <- if (missing(i) && missing(j)) {
                    x[[k]]
            }
            else if (is.matrix(x[[k]])) {
                if (missing(i)) {
                    x[[k]][, j, drop = drop]
                }
                else if (missing(j)) {
                    x[[k]][i, , drop = drop]
                }
                else {
                    x[[k]][i, j, drop = drop]
                }
            }
            else {
                ### Coerce x[[k]] to matrix before subsetting
                z <- as.matrix.mcmc(x[[k]])
                if (missing(i)) {
                    mcmc(z[, j, drop = TRUE], start(x), end(x), thin(x))
                }
                else if (missing(j)) {
                    z[i, , drop = TRUE]
                }
                else {
                    z[i, j, drop = TRUE]
                }
            }
        }
    }
    if (is.list(y) && all(sapply(y, is.mcmc, simplify = TRUE))) {
        y <- mcmc.list(y)
    }
    return(y)
}


"mcmc.list" <- function (...) 
{
  x <- list(...)
  if (length(x) == 1 && is.list(x[[1]])) 
    x <- x[[1]]
  if (!all(unlist(lapply(x, is.mcmc)))) 
    stop("Arguments must be mcmc objects")
  nargs <- length(x)
  if (nargs >= 2) {
    xmcpar <- lapply(x, mcpar)
    if (!all(unlist(lapply(xmcpar, "==", xmcpar[[1]])))) 
      stop("Different start, end or thin values in each chain")
    xnvar <- lapply(x, nvar)
    if (!all(unlist(lapply(xnvar, "==", xnvar[[1]])))) 
      stop("Different number of variables in each chain")
    xvarnames <- lapply(x, varnames, allow.null = FALSE)
    if (!all(unlist(lapply(xvarnames, "==", xvarnames[[1]])))) 
      stop("Different variable names in each chain")
  }
  if (is.R())
    class(x) <- "mcmc.list"
  else
    oldClass(x) <- "mcmc.list"
  return(x)
}

"start.mcmc.list" <- function (x, ...) 
{
  start(x[[1]])
}

"end.mcmc.list" <- function (x, ...) 
{
  end(x[[1]])
}

"thin.mcmc.list" <- function (x, ...) 
{
  thin(x[[1]])
}

"is.mcmc.list" <- function (x) 
  inherits(x, "mcmc.list")

"plot.mcmc.list" <-
  function (x, trace = TRUE, density = TRUE, smooth = TRUE, bwf, 
            auto.layout = TRUE, ask = par("ask"), ...) 
{
## RGA fixed to use default ask value.
  oldpar <- NULL
  on.exit(par(oldpar))
  if (auto.layout) {
    mfrow <- set.mfrow(Nchains = nchain(x), Nparms = nvar(x), 
                       nplots = trace + density)
    oldpar <- par(mfrow = mfrow)
  }
  for (i in 1:nvar(x)) {
    if (trace) 
      ## RGA fixed to propagate ... argument.
      traceplot(x[, i, drop = FALSE], smooth = smooth, ...)
    if (density) {
      if (missing(bwf)) 
        ## RGA fixed to propagate ... argument.
        densplot(x[, i, drop = FALSE], ...)
      else densplot(x[, i, drop = FALSE], bwf = bwf, ...)
    }
    if (i==1)
       oldpar <- c(oldpar, par(ask = ask))
  }
}


"summary.mcmc.list" <-
  function (object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) 
{
  x <- mcmc.list(object)
  statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
  varstats <- matrix(nrow = nvar(x), ncol = length(statnames), 
                     dimnames = list(varnames(x), statnames))
  xtsvar <- matrix(nrow = nchain(x), ncol = nvar(x))
  if (is.matrix(x[[1]])) {
    for (i in 1:nchain(x))
      for(j in 1:nvar(x))
        xtsvar[i, j] <- safespec0(x[[i]][,j])
    xlong <- do.call("rbind", x)
  }
  else {
    for (i in 1:nchain(x))
      xtsvar[i, ] <- safespec0(x[[i]])
    xlong <- as.matrix(x)
  }

  xmean <- apply(xlong, 2, mean)
  xvar <- apply(xlong, 2, var)
  xtsvar <- apply(xtsvar, 2, mean)
  varquant <- t(apply(xlong, 2, quantile, quantiles))
  varstats[, 1] <- xmean
  varstats[, 2] <- sqrt(xvar)
  ##RGA fixed so now give correct std error for pooled (across chains). 
  varstats[, 3] <- sqrt(xvar/(niter(x)*nchain(x)))
  varstats[, 4] <- sqrt(xtsvar/(niter(x)*nchain(x)))
  varquant <- drop(varquant)
  varstats <- drop(varstats)
  out <- list(statistics = varstats, quantiles = varquant, 
              start = start(x), end = end(x), thin = thin(x),
              nchain = nchain(x))
  class(out) <- "summary.mcmc"
  return(out)
}

"as.matrix.mcmc.list" <-
  function (x, iters = FALSE, chains = FALSE, ...) 
{
  x <- mcmc.list(x)
  y <- matrix(nrow = niter(x) * nchain(x), ncol = nvar(x) + 
              chains + iters)
  var.cols <- chains + iters + 1:nvar(x)
  for (i in 1:nchain(x)) {
    use.rows <- niter(x) * (i - 1) + 1:niter(x)
    if (chains) 
      y[use.rows, 1] <- i
    if (iters) 
      y[use.rows, chains + 1] <- as.vector(time(x))
    y[use.rows, var.cols] <- x[[i]]
  }
  rownames <- character(ncol(y))
  if (chains) 
    rownames[1] <- "CHAIN"
  if (iters) 
    rownames[1 + chains] <- "ITER"
  rownames[var.cols] <- varnames(x, allow.null = FALSE)
  dimnames(y) <- list(NULL, rownames)
  return(y)
}

"as.mcmc.mcmc.list" <- function (x, ...) 
{
  if (nchain(x) == 1) 
    return(x[[1]])
  else stop("Can't coerce mcmc.list to mcmc object:\n more than 1 chain")
}

"time.mcmc.list" <- function (x, ...) 
  time(x[[1]])

"window.mcmc.list" <- function (x, ...) 
{
  structure(lapply(x, window.mcmc, ...), class = "mcmc.list")
}

"head.mcmc.list" <- function (x, ...) 
{
  structure(lapply(x, head.mcmc, ...), class = "mcmc.list")
}


"tail.mcmc.list" <- function (x, ...) 
{
  structure(lapply(x, tail.mcmc, ...), class = "mcmc.list")
}

"as.mcmc.list" <- function (x, ...) 
  UseMethod("as.mcmc.list")

"as.mcmc.list.default" <- function (x, ...) 
  if (is.mcmc.list(x)) x else mcmc.list(x)

"as.array.mcmc.list" <- function(x, drop=TRUE, ...)
{
  y <- array(dim=c(niter(x), nvar(x), nchain(x)),
             dimnames = list(iter=time(x), var=varnames(x), chain=chanames(x)))
  for(i in 1:nchain(x))
    y[,,i] <- x[[i]]
  if(drop)
    return(drop(y))
  else
    return(y)
}

