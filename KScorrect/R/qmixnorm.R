#' @rdname dmixnorm
#' @export
#' @importFrom stats spline
qmixnorm <- function (p, mean, sd, pro, expand=1) {
  if(mode(p) != "numeric")
    stop("'p' must be a non-empty numeric vector")
  if (any(missing(mean), missing(sd)))
    stop("'mean' and 'sd' not provided, without default.")
  mean <- as.vector(mean, mode = "numeric")
  G <- length(mean)
  sd <- as.vector(sd, mode = "numeric")
  if (missing(pro)) {
    pro <- rep(1/G, G)
    warning("mixing proportion 'pro' not provided. Assigned equal proportions by default.")
  }
  if (any(pro < 0L, sd < 0L))
    stop("'pro' and 'sd' must not be negative.")
  lpro <- length(pro)
  lsd <- length(sd)
  if(lsd==1L & G > 1L) {
    sd[seq(G)] <- sd[1]
    lsd <- length(sd)
    warning("'equal variance model' implemented. If want 'variable-variance model', specify remaining 'sd's.")
  }
  if(G < lsd | G < lpro | (lsd > 1L & G != lsd) | (!missing(pro) & G != lpro))
    stop("the lengths of supplied parameters do not make sense.")
  pro <- as.vector(pro, mode = "numeric")
  pro <- pro/sum(pro)
  nr <- 10000
  x <- rmixnorm(nr * G, mean = mean, sd = sd, pro = pro)
  if(mode(expand) != "numeric" | expand < 0L)
    stop("'expand' must be a non-negative number.")
  span <- seq(min(x) - expand * diff(range(x)), max(x) + expand * diff(range(x)), length = nr)
  cdf <- vector(mode = "numeric", length = nr)
  for (g in seq.int(G)) {
    cdf <- cdf + pro[g] * pnorm(span, mean[g], sd[g])
  }
  quants <- stats::spline(cdf, span, method="fmm", xout=p)$y
  quants[which(p < 0L | p > 1L)] <- NaN
  quants[which(p == 0L)] <- -Inf
  quants[which(p == 1L)] <- Inf
  if(any(is.nan(quants)))
     warning("Some quantile values could not be calculated. If all 'p's are within [0,1], try reducing the value of 'expand' and try again.")
  return(as.vector(quants))
}
