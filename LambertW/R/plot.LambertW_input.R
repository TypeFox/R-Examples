#' @rdname LambertW_input_output-methods
#' @description
#' \code{plot.LambertW_input} plots the theoretical (1) pdf and (2) cdf of the
#' input \eqn{X \sim F_X(x \mid \boldsymbol \beta)}.
#' @export
plot.LambertW_input <- function(x, xlim = NULL, ...) {
  
  obj <- x
  tau <- obj$tau
  check_tau(tau)

  left.limit.is.zero <- FALSE
  if (obj$user.defined) {
    fam.type <- list(location = obj$tau["mu_x"] != 0,
                     scale = obj$tau["sigma_x"] != 1)
  } else {
    fam.type <- get_distname_family(obj$distname)
  }
  if (fam.type$location && fam.type$scale) {
    left.limit.is.zero <- FALSE
  } else if (!fam.type$location && fam.type$scale) {
    # if it's a scale family, but not location then the left 
    # limit should be zero (since its a non-negative)
    left.limit.is.zero <- TRUE
  }
  
  if (obj$user.defined) {
    coverage <- obj$q(c(0.005, 0.995))
  } else {
    coverage <- qLambertW(c(0.005, 0.995), theta = complete_theta(list(beta = obj$beta)), 
                          distname = obj$distname)
  }
  if (is.null(xlim)) {
    xlim <- coverage
    if (obj$distname == "chisq") {
      xlim[2] <- obj$tau["mu_x"] + (3 + 3 * (left.limit.is.zero)) * sqrt(2 * obj$beta)
    }
  }
  stopifnot(length(xlim) == 2,
            is.numeric(xlim),
            all(!is.na(xlim)),
            xlim[1] < xlim[2])
  
  x.seq <- seq(xlim[1], xlim[2], length = 101)
  ylim.pdf <- range(c(obj$d(x.seq)))  # y-range for density
  ylim.cdf <- range(c(obj$p(x.seq)))  # y-range for cdf
  
  if (obj$distname == "chisq" && obj$beta == 1) {
    ylim.pdf <- c(0, min(1, ylim.pdf[2]))
  }

  layout(matrix(1:2, nrow = 1))
  # pdf plot
  plot(obj$d, xlim[1], xlim[2], lwd = 1, main = obj$distname.with.beta, ylab = "pdf", 
       xlab = "x", ylim = ylim.pdf)
  grid()
  abline(v = obj$tau["mu_x"], lty = 2)
  # cdf plot
  plot(obj$p, xlim[1], xlim[2], lwd = 1, main = obj$distname.with.beta, ylab = "cdf", 
       xlab = "x", ylim = ylim.cdf)
  grid()
}