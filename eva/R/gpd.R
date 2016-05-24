#' The Generalized Pareto Distribution (GPD)
#'
#' Density, distribution function, quantile function and random number generation for the Generalized Pareto
#' distribution with location, scale, and shape parameters.
#' @name gpd
#' @rdname gpd
#' @param x Vector of observations.
#' @param q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param loc,scale,shape Location, scale, and shape parameters. Can be vectors, but
#' the lengths must be appropriate.
#' @param log.d Logical; if TRUE, the log density is returned.
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @examples
#' dgpd(2:4, 1, 0.5, 0.01)
#' dgpd(2, -2:1, 0.5, 0.01)
#' pgpd(2:4, 1, 0.5, 0.01)
#' qgpd(seq(0.9, 0.6, -0.1), 2, 0.5, 0.01)
#' rgpd(6, 1, 0.5, 0.01)
#'
#' ## Generate sample with linear trend in location parameter
#' rgpd(6, 1:6, 0.5, 0.01)
#'
#' ## Generate sample with linear trend in location and scale parameter
#' rgpd(6, 1:6, seq(0.5, 3, 0.5), 0.01)
#'
#' p <- (1:9)/10
#' pgpd(qgpd(p, 1, 2, 0.8), 1, 2, 0.8)
#' ## [1] 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
#'
#' ## Incorrect syntax (parameter vectors are of different lengths other than 1)
#' # rgpd(1, 1:8, 1:5, 0)
#'
#' ## Also incorrect syntax
#' # rgpd(10, 1:8, 1, 0.01)
#'
#' @details The Generalized Pareto distribution function is given (Pickands, 1975)
#' by \deqn{H(y) = 1 - \Big[1 + \frac{\xi (y - \mu)}{\sigma}\Big]^{-1/\xi}} defined
#' on \eqn{\{y : y > 0, (1 + \xi (y - \mu) / \sigma) > 0 \}}, with location \eqn{\mu},
#' scale \eqn{\sigma > 0}, and shape parameter \eqn{\xi}.
#'
#' @references Pickands III, J. (1975). Statistical inference using extreme order statistics. Annals of Statistics, 119-131.
NULL


#' @rdname gpd
#' @export
dgpd <- function(x, loc = 0, scale = 1, shape = 0, log.d = FALSE) {
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (length(x) > 1) &
    (((length(loc) != length(x)) & (length(loc) != 1)) |
    ((length(scale) != length(x)) & (length(scale) != 1)) |
    ((length(shape) != length(x)) & (length(shape) != 1)))
  cond2 <- (length(x) == 1) &
    (length(unique(c(length(x), length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  if(length(shape) == 1)
    shape <- rep(shape, max(length(x), length(loc), length(scale)))
  below.support <- x < loc
  x <- pmax(x, loc)
  x <- ifelse(shape >= 0, x, pmin(x, (loc - scale/shape)))
  w <- (x - loc) / scale
  log.density <- -log(scale) - ifelse(shape == 0, w, ((1/shape) + 1) * log1p(w * shape))
  log.density[is.nan(log.density) | is.infinite(log.density) | below.support] <- -Inf
  if(!log.d)
    log.density <- exp(log.density)
  log.density
}


#' @rdname gpd
#' @export
rgpd <- function(n, loc = 0, scale = 1, shape = 0) {
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (n > 1) &
    (((length(loc) != n) & (length(loc) != 1)) |
    ((length(scale) != n) & (length(scale) != 1)) |
    ((length(shape) != n) & (length(shape) != 1)))
  cond2 <- (n == 1) &
    (length(unique(c(n, length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  qgpd(runif(n), loc, scale, shape)
}


#' @rdname gpd
#' @export
qgpd <- function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  if(log.p)
    p <- exp(p)
  if((min(p, na.rm = TRUE) <= 0) || (max(p, na.rm = TRUE) >= 1))
    stop("`p' must contain probabilities in (0,1)")
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (length(p) > 1) &
    (((length(loc) != length(p)) & (length(loc) != 1)) |
    ((length(scale) != length(p)) & (length(scale) != 1)) |
    ((length(shape) != length(p)) & (length(shape) != 1)))
  cond2 <- (length(p) == 1) &
    (length(unique(c(length(p), length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  if(lower.tail)
    p <- 1 - p
  if(length(shape) == 1)
    shape <- rep(shape, max(length(p), length(loc), length(scale)))
  ifelse(shape == 0, loc - scale * log(p), loc + scale * expm1(-shape * log(p)) / shape)
}


#' @rdname gpd
#' @export
pgpd <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE, log.p = FALSE) {
  if(min(scale) <= 0)
    stop("Invalid scale")
  cond1 <- (length(q) > 1) &
    (((length(loc) != length(q)) & (length(loc) != 1)) |
    ((length(scale) != length(q)) & (length(scale) != 1)) |
    ((length(shape) != length(q)) & (length(shape) != 1)))
  cond2 <- (length(q) == 1) &
    (length(unique(c(length(q), length(loc), length(scale), length(shape)))) > 2)
  if(cond1 | cond2)
    stop("Invalid parameter length!")
  if(length(shape) == 1)
    shape <- rep(shape, max(length(q), length(loc), length(scale)))
  q <- pmax(q, loc)
  q <- ifelse(shape >= 0, q, pmin(q, (loc - scale/shape)))
  w <- (q - loc) / scale
  p <- ifelse(shape == 0, 1 - exp(-w), 1 - exp((-1/shape)*log1p(w*shape)))
  if(!lower.tail)
    p <- 1 - p
  if(log.p)
    p <- log(p)
  p
}

