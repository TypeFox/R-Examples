# Method for "density" generic:  Takes an spEM object and returns
# a corresponding KDE for the appropriate component and block, evaluated
# at the given set of points.  
# Does not use the FFT like the density.default function does; still
# quite fast, but not optimized for speed
density.spEM <- function (x, u = NULL, component=1, block = 1, 
                          scale = FALSE, ...) 
{  
  m <- NCOL(x$posteriors)
  r <- NCOL(x$data)
  n <- NROW(x$data)
  if (is.null(blockid <- x$blockid)) {
    coords <- 1
    u2 <- rep(1, r)
  } 
  else {
    u2 <- match(x$blockid, unique(x$blockid)) # convert blockid to integers 1, 2, ...
    coords <- blockid == block
    if (!any(coords)) 
      stop("Illegal value of block argument")
  }
  stackedx <- x$data[rep(1:n,m),]  
  cs <- colSums(x$posteriors)
  z.tmp <- sweep(x$posteriors, 2, cs, "/")
  z.tmp[, cs==0] <- 1/NROW(z.tmp) # Just in case
  wts <- rep(as.vector(z.tmp),r)
  mu <- matrix(x$muhat, nrow=m)
  sigma <- matrix(x$sigmahat, nrow=m)
  scaledx <- as.vector((stackedx - mu[rep(1:m, each=n), u2])/
                       sigma[rep(1:m, each=n), u2])
  bw <- x$bandwidth
  if (is.null(u)) {
    xx <- as.vector(as.matrix(x$data)[, coords])
    u = seq(min(xx) - 4 * bw, max(xx) + 4 * bw, len = 250)
  }
  # This part isn't used for now:
  if (!is.null(x$symmetric) && x$symmetric) {
    d <- wkde(x=scaledx, u=(u-mu[component, block])/sigma[component, block], 
              w=wts, bw=bw, sym=TRUE) / sigma[component, block]
  }
  else {
    d <- wkde(x=scaledx, u=(u-mu[component, block])/sigma[component, block], 
              w=wts, bw=bw) / sigma[component, block]
  }
  if (scale) 
    d <- d * x$lambdahat[component]
  structure(list(x = u, y = d, bw = bw, n = n, call = match.call(), 
                 data.name = deparse(substitute(x)), has.na = FALSE), 
             class = "density")
}

