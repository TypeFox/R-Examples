kde.quantiles <- function(x, probs = c(0.05, 0.95), ...) 
{
  if(is.null(x))
    stop("nothing to do!")
  if(is.null(probs))
    probs <- c(0.05, 0.95)
  args <- list(...)
  if(is.null(args$from))
    args$from <- min(x, na.rm = TRUE)
  if(is.null(args$to))
    args$to <- max(x, na.rm = TRUE)
  args$x <- x
  kde <- do.call(stats::density.default, delete.args(stats::density.default, args))
  quants <- approx.quantile(x = kde$x, y = kde$y, probs = probs)

  return(quants)
}


approx.quantile <- function(x, y, probs = c(0.05, 0.95))
{
  afun <- stats::splinefun(x = x, y = y)
  foo <- function(x, afun, minx, prob) {
    if(x < minx)
      x <- minx + 1
    quant <- stats::integrate(f = afun, lower = minx, upper = x)$value
    return((quant - prob)^2)
  }
  minx = min(x)
  rangex <- range(x)
  n <- length(probs)
  rval <- rep(NA, length = n)
  for(k in 1:n) {
    rval[k] <- stats::optimize(f = foo, interval = rangex, afun = afun, 
      minx = minx, prob = probs[k])$minimum
  }

  return(rval)
}

