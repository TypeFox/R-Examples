#' @title Discovery Curve
#' @description Calculate the components of a species discovery curve.
#' 
#' @param f a vector of species frequencies where \code{f[i]} is the number of 
#'   species represented by only \code{i} samples.
#' @param max.x the maximum number of samples to calculate the curve for. 
#'   Defaults to the sample size of \code{f}.
#' @param n.pts number of points between 0 and \code{max.x} to estimate.
#' @param ci size of the confidence interval (0.5:1).
#' @param f0.func function to use to calculate \code{\link{f0}}.
#' @param ... other arguments to \code{f0.func}.
#' 
#' @return a list with:
#' \item{f.stats}{a named vector from \code{f0.func}.}
#' \item{s.ind}{a \code{matrix} of S.ind estimates for each value of m along 
#'   with the standard deviation of S.ind.}
#' \item{s.ind.ci}{a \code{matrix} of the upper and lower confidence intervals 
#'   of S.ind.}
#' \item{ci.poly}{a \code{matrix} of points describing the ci polygon.}
#' \item{rarefact.line}{a \code{matrix} of points defining the rarefaction 
#'   line (<= S.obs).}
#' \item{extrap.line}{a \code{matrix} of points defining the extrapolation 
#'   line (> S.obs).}
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references Colwell, R.K., A. Chao, N.J. Gotelli, S.-Y. Lin, C.X. Mao, 
#'   R.L. Chazdon, and J.T. Longino. 2012. Models and estimators linking 
#'   individual-based and sample-based rarefaction, extrapolation and 
#'   comparison of assemblages. Journal of Plant Ecology 5(1):3-21.
#' 
#' @seealso \code{\link{plot.discovery.curve}}
#' 
#' @examples
#' data(osa.old.growth)
#' f <- expand.freqs(osa.old.growth)
#' d <- discovery.curve(f, f0.func = Chao1, max.x = 1200)
#' plot(d)
#' 
#' @importFrom stats qnorm
#' @export
#' 
discovery.curve <- function(f, f0.func, max.x = sum(f * 1:length(f)), 
                            n.pts = 100, ci = 0.95, ...) {
  n <- sum(f * 1:length(f))
  n.seq <- ceiling(seq(0, ceiling(max.x), length.out = n.pts))
  n.seq <- sort(unique(c(n, n.seq)))
  s.ind <- t(sapply(n.seq, expected.num.species, f = f, f0.func = f0.func, ...))
  f.stats <- s.ind[1, 4:7]
  s.ind <- s.ind[, -(4:7)]
  s.ind <- s.ind[, c(3, 1, 2)]
  
  s.ind.ci <- s.ind[, "s.ind"] + qnorm((1 - ci) / 2) * 
    cbind(s.ind[, "sd.s.ind"], -s.ind[, "sd.s.ind"])
  ci.poly <- cbind(x = c(n.seq, rev(n.seq)), 
                   y = c(s.ind.ci[, 1], rev(s.ind.ci[, 2])))
  ci.poly[, "y"] <- ifelse(ci.poly[, "y"] < 1, 1, ci.poly[, "y"])
  i <- which(s.ind[, "m"] <= f.stats["n"])
  rarefact.line <- cbind(x = s.ind[i, "m"], y = s.ind[i, "s.ind"])
  extrap.line <- cbind(x = s.ind[-i, "m"], y = s.ind[-i, "s.ind"])
  
  d <- list(f.stats = f.stats, s.ind = s.ind, s.ind.ci = s.ind.ci, 
            ci.poly = ci.poly, rarefact.line = rarefact.line, 
            extrap.line = extrap.line)
  class(d) <- c("discovery.curve", class(d))
  return(d)
}
