#' @rdname LambertW_fit-methods
#' 
#' @description
#' \code{plot.LambertW_fit} plots a (1) histogram, (2) empirical density of the
#' data \code{y}. These are compared (3) to the theoretical \eqn{F_X(x \mid
#' \widehat{\boldsymbol \beta})} and (4) Lambert W \eqn{\times} 
#' \eqn{F_X(y \mid \widehat{\boldsymbol \beta})} densities.
#' 
#' @param xlim lower and upper limit of x-axis for cdf and pdf plots.
#' @param show.qqplot should a Lambert W\eqn{ \times} F QQ plot be displayed? Default: \code{FALSE}.
#' @export
#' 
plot.LambertW_fit <- function(x, xlim = NULL, show.qqplot = FALSE, ...) {
  
  xlim.by.user <- xlim
  yy <- x$data
  tau.complete <- complete_tau(x$tau)

  if (identical(x$method, "IGMM")) {
    x$distname <- "normal"
    x$theta <- tau2theta(x$tau, beta = x$tau[c("mu_x", "sigma_x")])
  } else if (identical(x$distname, "normal")) {
    x$beta <- x$tau[c("mu_x", "sigma_x")]
  }
  
  beta.y <- estimate_beta(yy, distname = x$distname)
  
  .PdfLambertW <- function(xx) {
    return(dLambertW(xx, theta = x$theta, distname = x$distname,
                     use.mean.variance = x$use.mean.variance))
  }
  .PdfZeroLambertW = function(xx){
    return(dLambertW(xx, theta = list(beta = beta.y), 
                     distname = x$distname,
                     use.mean.variance = x$use.mean.variance))
  }
  
  coverage <- qLambertW(c(0.005, 0.995), 
                        theta = x$theta, 
                        distname = x$distname,
                        use.mean.variance = x$use.mean.variance)
  x.lower <- coverage[1]
  x.upper <- coverage[2]
  if (!is.null(xlim.by.user[1])) {
    x.lower <- xlim.by.user[1]
  }
  if (!is.null(xlim.by.user[2])) {
    x.upper <- xlim.by.user[2]
  }
  
  COL <- c(1, 2, 4)  # Kernel, LambertW, Original
  LWD <- c(2, 3, 1)
  LTY <- c(1, 1, 2)

  # get good breaks for histogram
  good.num.breaks <- .optimalNumberOfBinsForHist(yy)

  dist.family <- get_distname_family(x$distname)
  if (dist.family$is.non.negative) {
    pdf.kde <- density(yy, from = 0)
  } else {
    pdf.kde <- density(yy)
  }
  hist.est <- hist(yy, good.num.breaks, plot = FALSE)
  pdf.para <- try(.PdfZeroLambertW(seq(x.lower, x.upper, length = 100)),
                  silent = TRUE)
  if (inherits(pdf.para, "try-error")) {
    pdf.para <- rep(0, 100)
  }

  # find support
  rv.support <- get_support(x$tau)
  sup.l <- seq(max(rv.support[1], x.lower), min(rv.support[2], x.upper), length = 100)
  pdf.paraLW <- .PdfLambertW(sup.l[-c(1, length(sup.l))])
  
  leg.txt <- c("Kernel", 
               paste0("Lambert W x ", x$distname, "\n (type: '", x$type, "')"), 
               x$distname)
  
  if (x$type == "s") {
    skewness.pos <- x$theta$gamma > 0
  } else if (x$type == "hh") {
    skewness.pos <- x$theta$delta[1] < x$theta$delta[2]
  } else {
    skewness.pos <- skewness(x$data) > 0
  }
  legend.pos <- ifelse(skewness.pos, "topright", "topleft")
  
  y.lim <- range(pdf.kde$y, pdf.para, pdf.paraLW, hist.est$intensities) * 1.1
  hist(yy, good.num.breaks, xlim = c(x.lower, x.upper),
       ylim = y.lim, prob = TRUE, density = 20, col = "darkgray", 
       ylab = "", main = "", ...)
  legend(legend.pos, leg.txt, col = COL, lwd = LWD, lty = LTY, cex = 0.8,
         box.lty = 0, horiz = FALSE, adj = 0)
  
  lines(pdf.kde, lwd = LWD[1], lty = LTY[1])
  if (any(pdf.para > 0)) {
    plot(.PdfZeroLambertW, x.lower, x.upper, add = TRUE, 
         lty = LTY[3], col = COL[3], lwd = LWD[3])
  }
  plot(.PdfLambertW, min(sup.l), max(sup.l), 
       lwd = LWD[2], col = COL[2], lty = LTY[2], ylab = "", add = TRUE)
  abline(v = rv.support, lwd = 2, lty = 3, col = COL[2])
  
  dist.family <- get_distname_family(x$distname)
  if (x$type == "s" && !dist.family$is.non.negative) {
    if (x$tau["gamma"] > 0) {
      legend.pos <- "right"
      bound.at <- rv.support[1]
      bound.legend.title <- "Lower bound"
    } else {
      legend.pos <- "left"
      bound.at <- rv.support[2]
      bound.legend.title <- "Upper bound"
    }
    par.tmp <- par()
    bound.legend.txt <- paste(names(bound.at), "=", round(bound.at, 2))
    
    text(ifelse(x$tau["gamma"] < 0, x.upper, x.lower), par.tmp$yaxp[2] * 0.9, 
         paste0(bound.legend.title, "\n", bound.legend.txt),
         pos = ifelse(x$tau["gamma"] < 0, 2, 4),
         cex = 0.75)
    #mtext(paste0(bound.legend.title, "\n", bound.legend.txt), side = 1, line = 0)
    #legend(legend.pos, bound.legend.txt, title = bound.legend.title, box.lty = 0)
  }
  box(lwd = 2)
  
  if (show.qqplot) {
    qqLambertW(yy, theta = x$theta, distname = x$distname,
               use.mean.variance = x$use.mean.variance)
  }
}