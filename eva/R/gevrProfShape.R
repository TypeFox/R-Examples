#' GEVr Shape Parameter Profile Likelihood Estimation for Stationary Models
#'
#' Computes the profile likelihood based confidence interval for the shape parameter of the stationary GEVr model.
#'
#' @param z A class object returned from gevrFit.
#' @param conf Confidence level to use. Defaults to 95 percent.
#' @param opt Optimization method to maximize the profile likelihood, passed to optim. The default method is Nelder-Mead.
#'
#' @examples
#' ## Compare the length of the shape confidence intervals using GEV1 vs. GEV10
#' set.seed(7)
#' x <- rgevr(200, 10, loc = 0.5, scale = 1, shape = 0.3)
#' z1 <- gevrFit(x[, 1])
#' z2 <- gevrFit(x)
#' gevrProfShape(z1)
#' gevrProfShape(z2)
#' @return
#' \item{Estimate}{Estimated shape parameter.}
#' \item{CI}{Profile likelihood based confidence interval for the shape parameter.}
#' \item{ConfLevel}{The confidence level used.}
#' @export
gevrProfShape <- function(z, conf = .95, opt = c("Nelder-Mead")) {
  if(z$gumbel | !z$stationary)
    stop("Object cannot be from a Gumbel and/or a nonstationary fit!")
  data <- as.matrix(z$data)
  theta <- as.numeric(z$par.ests)
  opt <- match.arg(opt)
  sol <- c(theta[1], theta[2])
  gevrLikShape <- function(a, sh) {
    if(a[2] <= 0) {
      out <- .Machine$double.xmax
    } else {
      out <- dgevr(data, loc = a[1], scale = a[2], shape = sh, log.d = TRUE)
      out <- - sum(out)
      if(out == Inf)
        out <- .Machine$double.xmax
    }
    out
  }
  cutoff <- qchisq(conf, 1)
  prof <- function(sh) {
    lmax <- dgevr(data, theta[1], theta[2], theta[3], log.d = TRUE)
    lmax <- sum(lmax)
    yes <- optim(sol, gevrLikShape, method = opt, sh = sh)
    sol <- yes$par
    lci <- -yes$value
    2*(lmax-lci) - cutoff
  }
  prof <- Vectorize(prof)
  suppressWarnings(out1 <- uniroot(prof, c(theta[3] - 1e-6, theta[3]), extendInt="downX"))
  suppressWarnings(out2 <- uniroot(prof, c(theta[3], theta[3] + 1e-6), extendInt="upX"))
  CI <- c(min(out1$root, out2$root), max(out1$root, out2$root))
  out <- list(theta[3], CI, conf)
  names(out) <- c("Estimate", "CI", "ConfLevel")
  out
}



