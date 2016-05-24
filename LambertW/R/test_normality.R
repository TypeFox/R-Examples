#' @title Visual and statistical Gaussianity check
#' 
#' @description
#' Graphical and statistical check if data is Gaussian (three common Normality
#' tests, QQ-plots, histograms, etc).
#' 
#' \code{test_normality} does not show the autocorrelation function (ACF)
#' estimate for lag \eqn{0}, since it always equals \eqn{1}. Thus removing it
#' does not lose any information, but greatly improves the y-axis scale for
#' higher order lags (which are usually very small compared to 1).
#' 
#' @param data a numeric vector of data values.
#' @param plot Should visual checks (histogram, densities, qqplot, ACF) be
#'     plotted?  Default \code{TRUE}; otherwise only hypothesis test results are
#'     returned.
#' @param show.volatility logical; if \code{TRUE} the squared (centered) data
#'     and its ACF are also shown. Useful for time series data to see if squares
#'     exhibit dependence (for financial data they typically do); default:
#'     \code{FALSE}.
#' @param pch a vector of plotting characters or symbols; default \code{pch =
#'     1}.
#' @param add.legend logical; if \code{TRUE} (default) a legend is placed in
#'     histogram/density plot.
#' @param seed optional; if sample size > 5,000, then some normality tests fail
#'     to run.  In this case it uses a subsample of size 5,000.  For
#'     reproducibility, the seed can be specified by user.  By default it uses a
#'     random seed.
#' @return A list with results of 3 normality tests (each of class \code{htest})
#'     and the seed used for subsampling: \item{anderson.darling}{Anderson
#'     Darling (if \pkg{nortest} package is available),}
#'     \item{shapiro.francia}{Shapiro-Francia (if \pkg{nortest} package is
#'     available),} \item{shapiro.wilk}{Shapiro-Wilk,} \item{seed}{seed for
#'     subsampling (only used if sample size > 5,000).}
#' 
#' @seealso \code{\link[stats]{shapiro.test}} in the \pkg{stats} package;
#'     \code{\link[nortest]{ad.test}}, \code{\link[nortest]{sf.test}} in the
#'     \pkg{nortest} package.
#' @references Thode Jr., H.C. (2002): \dQuote{Testing for Normality}. Marcel
#'     Dekker, New York.
#' @keywords htest hplot
#' @export
#' @examples
#' 
#' y <- rLambertW(n = 1000, theta = list(beta = c(3, 4), gamma = 0.3),
#'                distname = "normal")
#' test_normality(y)
#' 
#' x <- rnorm(n = 1000)
#' test_normality(x)
#' 
#' # mixture of exponential and normal
#' test_normality(c(rexp(100), rnorm(100, mean = -5)))
#' 
test_normality <- function(data, show.volatility = FALSE, plot = TRUE, pch = 1, 
                    add.legend = TRUE, seed = sample(1e6, 1)) {
  if (!is.numeric(data)) {
    data <- data$res
  }
  
  stopifnot(is.numeric(unlist(data)))
  
  data.label <- deparse(substitute(data))
  num.samples <- length(data)
  
  ## show visual checks for Normality
  if (plot) {
    if (show.volatility) {
      layout.dim <- c(3, 2)
      data <- ts(data)
    } else {
      layout.dim <- c(2, 2)
    }
    
    par(mfrow = layout.dim, mar = c(4.5, 4.5, 2, 0.5))
    plot(data, pch = pch, ylab = data.label)
    grid()
    
    #########
    .hist_and_density <- function(yy) {
      .PdfNorm <- function(xx) {
        dnorm(xx, mean = mu.y, sd = sd.y)
      }
      
      yy.range <- range(yy)
      x.lower <- yy.range[1] - 0.05 * abs(yy.range[1])
      x.upper <- yy.range[2] + 0.05 * abs(yy.range[2])
      
      mu.y <- mean(yy)
      sd.y <- sd(yy)
      
      COL <- c("kernel" = 1, "normal" =  2)
      LWD <- c("kernel" = 2, "normal" = 2)
      LTY <- c("kernel" = 1, "normal" = 2)
      num.breaks <- .optimalNumberOfBinsForHist(yy)
      hist.est <- hist(x = yy, breaks = num.breaks, plot = FALSE)
      pdf.kde <- density(yy)$y
      pdf.para <- .PdfNorm(seq(x.lower, x.upper, length = 100))
      
      hist(yy, breaks = num.breaks, 
           xlim = c(x.lower, x.upper), 
           ylim = range(pdf.kde, pdf.para, hist.est$intensities) * 1.2, 
           prob = TRUE, density = 10, col = "darkgray", 
           main = "Density estimates", 
           ylab = "", xlab = data.label)
      lines(density(yy), 
            lwd = LWD["kernel"], col = COL["kernel"], lty = LTY["kernel"])
      plot(.PdfNorm, x.lower, x.upper, add = TRUE, 
           lty = LTY["normal"], col = COL["normal"], lwd = LWD["normal"])
      
      if (add.legend) {
        legend.pos <- ifelse(skewness(yy) >= 0, "topright", "topleft")
        legend2.pos <- ifelse(skewness(yy) >= 0, "topleft", "topright")
        
        legend(legend.pos, 
               c(paste("skewness:", format(skewness(yy), digits = 2)), 
                 paste("kurtosis:  ", format(kurtosis(yy), digits = 2))),
               box.lty = 0)
        legend(legend2.pos,
               legend = sapply(c(bquote(N(.(round(mu.y, 2)), .(round(sd.y, 2))^2)),
                                 "Kernel"), as.expression),
               lty = LTY[c("normal", "kernel")], 
               col = COL[c("normal", "kernel")], 
               lwd = LWD[c("normal", "kernel")], box.lty = 0)
      }
    }
    ########
    .acf_no_zero <- function(x, ...) {
      acf.tmp <- acf(x, plot = FALSE, ...)
      CI.vals <- c(-1, 1) * qnorm(1 - 0.05/2) / sqrt(length(x))
      acf(x, ylim = range(acf.tmp$acf[-1], CI.vals) * 1.1, 
          xlim = c(1.25, length(acf.tmp$acf) - 1), ...)
      abline(h = CI.vals, lty = 2, col = "blue")
    }
    
    .acf_no_zero(data)
    .hist_and_density(data)
    qqnorm(data)
    qqline(data)
    grid()
    
    if (show.volatility) {
      data.centered <- scale(data, center = TRUE, scale = FALSE)
      if (is.ts(data)) {
        data.centered <- ts(data.centered)
      }
      plot(data.centered^2, ylab = paste0("Squared (centered) \n", 
                                          data.label))
      grid()
      .acf_no_zero(data.centered^2)
    }
  }
  
  # initialize output with NA
  out <- list(seed = seed,
              shapiro.wilk = NA,
              shapiro.francia = NA,
              anderson.darling = NA)

  if (requireNamespace("nortest", quietly = TRUE)) {
    # Run the normality hypothesis tests; 
    # Anderson-Darling does not have a restriction on sample size.
    out[["anderson.darling"]] <- nortest::ad.test(data)
  }
  # only update tests if 'nortest' package is available; otherwise return 'NA'
  if (num.samples > 5000) {
    cat("Shaprio-Wilk and Shapiro-Francia tests cannot be computed for ",
        "sample size > 5000.\n", 
        "Use random subsample of size 5000 instead.\n")
    set.seed(seed)
    data.test <- sample(data, size = 5000)
  } else {
    data.test <- data
  }
  if (requireNamespace("nortest", quietly = TRUE)) {
    out[["shapiro.francia"]] <- nortest::sf.test(data.test)
  }
  out[["shapiro.wilk"]] <- shapiro.test(data.test)
  
  return(out)
} 

#' @rdname test_normality
#' @export
#' @param ... arguments as in \code{test_normality}.
#' @description
#' \code{test_norm} is a shortcut for \code{test_normality}.

test_norm <- function(...) {
  test_normality(...)
}
