#' GEVr Return Level Estimate and Confidence Interval for Stationary Models
#'
#' Computes stationary m-period return level estimate and interval, using either the delta method or profile likelihood.
#'
#' @param z A class object returned from gevrFit. Must be a stationary fit.
#' @param period The number of periods to use for the return level.
#' @param conf Confidence level. Defaults to 95 percent.
#' @param method The method to compute the confidence interval - either delta method (default) or profile likelihood.
#' @param opt Optimization method to maximize the profile likelihood if that is selected. The default method is Nelder-Mead.
#'
#' @details It is generally accepted that profile likelihood confidence intervals provide greater accuracy than the delta
#' method, in particular for large return level periods. Also, by their nature, delta method confidence intervals must be symmetric
#' which may be undesirable for return level estimation. If the original fit was Gumbel, then return levels will be for the Gumbel
#' distribution.
#' @references http://www.mas.ncl.ac.uk/~nlf8/teaching/mas8391/background/chapter2.pdf
#' @references Coles, S. (2001). An introduction to statistical modeling of extreme values (Vol. 208). London: Springer.
#' @examples
#' x <- rgevr(100, 2, loc = 0.5, scale = 1, shape = 0.3)
#' z <- gevrFit(x)
#' ## Compute 250-period return level.
#' gevrRl(z, 250, method = "delta")
#' @return
#' \item{Estimate}{Estimated m-period return level.}
#' \item{CI}{Confidence interval for the m-period return level.}
#' \item{Period}{The period length used.}
#' \item{ConfLevel}{The confidence level used.}
#' @details Caution: The profile likelihood optimization may be slow (on the order of minutes).
#' @export
gevrRl <- function(z, period, conf = .95, method = c("delta", "profile"), opt = c("Nelder-Mead")) {
  if(!z$stationary)
    stop("Return levels can only be produced for the stationary model!")
  method <- match.arg(method)
  data <- as.matrix(z$data)
  theta <- z$par.ests
  cov <- z$varcov
  m <- -log1p(-(1/period))
  if(!z$gumbel) est <- theta[1] - (theta[2]/theta[3])*(-expm1(-theta[3]*log(m))) else est <- theta[1] - theta[2]*log(m)
  est <- as.numeric(est)
  if(method == "delta") {
    if(!z$gumbel) {
      del <- matrix(ncol = 1, nrow = 3)
      del[1, 1] <- 1
      del[2, 1] <- -((theta[3])^(-1))*(-expm1(-theta[3]*log(m)))
      del[3, 1] <- ((theta[2])*(theta[3]^(-2))*(-expm1(-theta[3]*log(m)))) - ((theta[2])*((theta[3])^(-1))*((1+expm1(-theta[3]*log(m)))*log(m)))
      se <- as.numeric(sqrt(t(del) %*% cov %*% del))
    } else {
      se <- as.numeric(sqrt(cov[1, 1] - ((cov[2, 2] + cov[2, 1]) * log(m)) + (cov[2, 2] * (log(m))^2)))
    }
    alpha <- (1-conf)/2
    lower <- est - qnorm(1-alpha)*se
    upper <- est + qnorm(1-alpha)*se
    CI <- c(lower, upper)
  } else {
    opt <- match.arg(opt)
    if(!z$gumbel) sol <- c(theta[2], theta[3]) else sol <- c(theta[2])
    gevrLik <- function(a, xp) {
      loc <- xp + (a[1]/a[2])*(-expm1(-a[2]*log(m)))
      if(a[1] <= 0) {
       out <- .Machine$double.xmax
      } else {
        out <- dgevr(data, loc = loc, scale = a[1], shape = a[2], log.d = TRUE)
        out <- - sum(out)
        if(out == Inf)
          out <- .Machine$double.xmax
      }
      out
    }
    gumLik <- function(a, xp) {
      loc <- xp + a*log(m)
      if(a <= 0) {
        out <- .Machine$double.xmax
      } else {
        out <- dgevr(data, loc = loc, scale = a, shape = 0, log.d = TRUE)
        out <- - sum(out)
        if(out == Inf)
          out <- .Machine$double.xmax
      }
      out
    }
    cutoff <- qchisq(conf, 1)
    prof <- function(xp) {
      if(!z$gumbel) {
        lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
        lmax <- sum(lmax)
        yes <- optim(sol, gevrLik, method = opt, xp = xp)
      } else {
        lmax <- dgevr(data, theta[1], theta[2], 0, log.d = TRUE)
        lmax <- sum(lmax)
        yes <- optim(sol, gumLik, method = opt, xp = xp)
      }
      sol <- yes$par
      lci <- -yes$value
      2*(lmax-lci) - cutoff
    }
    prof <- Vectorize(prof)
    suppressWarnings(out1 <- uniroot(prof, c(est - 1e-6, est), extendInt="downX"))
    suppressWarnings(out2 <- uniroot(prof, c(est, est + 1e-6), extendInt="upX"))
    CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  }
  out <- list(est, CI, period, conf)
  names(out) <- c("Estimate", "CI", "Period", "ConfLevel")
  out
}
