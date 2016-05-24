# Method for "density" generic:  Takes an npEM object and returns
# a corresponding KDE for the appropriate component and block, evaluated
# at the given set of points.  
# Does not use the FFT like the density.default function does; still
# quite fast, but not optimized for speed
density.npEM <- function (x, u = NULL, component = 1, block = 1, scale = FALSE, 
    ...) 
{
  if (is.null(blockid <- x$blockid)) {
    coords <- 1
  } 
  else {
    coords <- blockid == block
    if (!any(coords)) 
      stop("Illegal value of block argument")
  }
  m <- NCOL(x$post)
  xx <- as.vector(as.matrix(x$data)[, coords])
  if(is.matrix(x$bandwidth))
    bw <- x$bandwidth[block, component]
  else bw <- x$bandwidth
  if (is.null(u)) {
    u = seq(min(xx) - 4 * bw, max(xx) + 4 * bw, len = 250)
  }
  if (component > m || component < 1) 
    stop("Illegal value of component argument")
  if (!is.null(x$symmetric) && x$symmetric) {
    n <- length(xx)
    d <- wkde(x = rep(xx, m) - rep(x$muhat, each = n), u = u - 
              x$muhat[component], w = as.vector(x$post), bw = bw, 
              sym = TRUE)
  }
  else {
    n <- NROW(x$data)
    wts <- rep(x$post[, component], sum(coords))
    d <- wkde(xx, u = u, w = wts, bw = bw)
  }
  if (scale) 
    d <- d * x$lambdahat[component]
  structure(list(x = u, y = d, bw = bw, n = n, call = match.call(), 
                 data.name = deparse(substitute(x)), has.na = FALSE), 
             class = "density")
}

