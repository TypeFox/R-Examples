#' @rdname LambertW_input_output-methods
#' @description
#' \code{plot.LambertW_output} plots the theoretical (1) pdf and (2) cdf of the
#' output RV \eqn{Y \sim} Lambert W \eqn{\times} \eqn{F_X(x \mid \boldsymbol
#' \beta)}. It overlays the plot with the pdf and cdf of the input RV \eqn{X \sim F_X(x \mid
#' \boldsymbol \beta)} (setting \eqn{\gamma = \delta = 0, \alpha = 1}).
#' @export

plot.LambertW_output <- function(x, xlim = NULL, ...) {
  
  obj <- x
  check_theta(obj$theta, distname = x$input.distname)
  check_tau(obj$tau)
  
  if (is.null(xlim)) {
    xlim <- qLambertW(c(0.005, 0.995), theta = x$theta, distname = x$input.distname)
  }
  left.limit.is.zero <- get_distname_family(obj$input.distname)$scale && 
    !get_distname_family(obj$input.distname)$location
  
  # set limit to zero if its a non-negative distribution, and if the density does not
  # go to infty at 0
  if (left.limit.is.zero && obj$d(max(xlim[1] / 10, 1e-5)) < 100) {
    xlim[1] <- 0
  }
  
  stopifnot(length(xlim) == 2,
            is.numeric(xlim),
            all(!is.na(xlim)),
            xlim[1] < xlim[2])
  
  x.seq <- seq(xlim[1], xlim[2], length = 101)
  
  input.distname.with.beta <- gsub("Lambert W x ", "", obj$distname.with.beta)
  input.distname <- gsub("Lambert W x ", "", obj$distname)

  # params for input distribution (set gamma = delta = 0 and alpha = 1 by default)
  params.input <- complete_theta(list(beta = obj$theta$beta))
  
  ylim.pdf <- range(c(0, obj$d(x.seq, obj$theta), obj$d(x.seq, params.input)))  # y-range for density
  ylim.cdf <- range(c(obj$p(x.seq, obj$theta), obj$p(x.seq, params.input)))  # y-range for cdf
  
  if (obj$input.distname == "chisq" && obj$theta$beta == 1) {
    ylim.pdf <- c(0, min(1, ylim.pdf[2]))
  }
  
  layout(matrix(1:2, nrow = 1))
  # pdf plot
  curve(obj$d(x, obj$theta), xlim[1], xlim[2], lwd = 2, lty = 2, col = 2, 
        main = obj$distname.with.beta, 
        ylab = "pdf", xlab = "y", 
        ylim = c(max(0, ylim.pdf[1]), ylim.pdf[2] * 1.25))
  grid()
  curve(obj$d(x, params.input), xlim[1], xlim[2], add = TRUE, lty = 1, col = 1, lwd = 2)
  legend("topright", lwd = 2, col = 2:1, lty = 2:1, cex = 1, 
         c(obj$distname.with.beta, input.distname.with.beta),
         box.lty = 0)

  if (obj$type %in% c("s", "h")) {
    mtext(substitute(list(alpha == al, gamma == b, delta == a), 
                     list(a = round(obj$theta$delta, 3), 
                          al = round(obj$theta$alpha, 3), b = round(obj$theta$gamma, 3))))
  } else if (x$type == "hh") {
    mtext(substitute(list(alpha == al, gamma = g, delta[l] == a, delta[r] == b), 
                     list(g = round(obj$theta$gamma, 3), al = round(obj$theta$alpha, 3), 
                          a = round(obj$theta$delta, 3)[1], b = round(obj$theta$delta, 3)[2])))
  }
  
  # cdf plot
  curve(obj$p(x, obj$theta), xlim[1], xlim[2], lwd = 2, lty = 2, col = 2, 
        main = obj$distname.with.beta, 
        ylab = "cdf", xlab = "y", 
        ylim = c(max(0, ylim.cdf[1]), ylim.cdf[2] * 1.25))
  grid()
  curve(obj$p(x, params.input), xlim[1], xlim[2], add = TRUE, lty = 1, col = 1, lwd = 2)
  abline(h = c(0, 1))
  legend("topleft", lwd = 2, col = 2:1, lty = 2:1, cex = 1, 
         c(obj$distname.with.beta, input.distname.with.beta),
         box.lty = 0)

  if (x$type %in% c("s", "h")) {
    mtext(substitute(list(alpha == al, gamma == b, delta == a), 
                     list(a = round(obj$theta$delta, 3), al = round(obj$theta$alpha, 3), 
                          b = round(obj$theta$gamma, 3))))
  } else if (x$type == "hh") {
    mtext(substitute(list(alpha == al, gamma = g, delta[l] == a, delta[r] == b), 
                     list(g = round(obj$theta$gamma, 3), al = round(obj$theta$alpha, 3), 
                          a = round(obj$theta$delta, 3)[1], b = round(obj$theta$delta, 3)[2])))
  }
}
  