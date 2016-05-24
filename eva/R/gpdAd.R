gpdAdGen <- function(n, theta) {
  data1 <- rgpd(n, loc = 0, scale = theta[1], shape = theta[2])
  fit1 <- tryCatch(gpdFit(data1, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit1)) {
    teststat <- NA
  } else {
    scale1 <- fit1$par.ests[1]
    shape1 <- fit1$par.ests[2]
    thresh1 <- findthresh(data1, n)
    newdata1 <- pgpd(data1, loc = thresh1, scale = scale1, shape = shape1)
    newdata1 <- sort(newdata1)
    i <- seq(1, n, 1)
    teststat <- -n - (1/n)*sum((2*i - 1)*(log(newdata1) + log1p(-rev(newdata1))))
  }
  teststat
}


#' Generalized Pareto Distribution Anderson-Darling Test
#'
#' Anderson-Darling goodness-of-fit test for the Generalized Pareto (GPD) distribution.
#' @param data Data should be in vector form, assumed to be from the GPD.
#' @param bootstrap Should bootstrap be used to obtain p-values for the test? By default, a table of critical values is used via interpolation. See details.
#' @param bootnum Number of replicates if bootstrap is used.
#' @param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @references Choulakian, V., & Stephens, M. A. (2001). Goodness-of-fit tests for the Generalized Pareto distribution. Technometrics, 43(4), 478-484.
#' @examples
#' ## Generate some data from GPD
#' x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
#' gpdAd(x)
#' @details A table of critical values were generated via Monte Carlo simulation for shape
#' parameters -0.5 to 1.0 by 0.1, which provides p-values via log-linear interpolation from
#' .001 to .999. For p-values below .001, a linear equation exists by regressing -log(p-value)
#' on the critical values for the tail of the distribution (.950 to .999 upper percentiles). This
#' regression provides a method to extrapolate to arbitrarily small p-values.
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Estimated value of theta for the initial data.}
#' \item{effective_bootnum}{Effective number of bootstrap replicates if bootstrap
#' based p-value is used (only those that converged are used).}
#' @import parallel
#' @export

gpdAd <- function(data, bootstrap = FALSE, bootnum = NULL, allowParallel = FALSE, numCores = 1) {
  if(bootstrap == TRUE & is.null(bootnum))
    stop("Must specify some number of boostrap samples")
  n <- length(data)
  fit <- tryCatch(gpdFit(data, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit))
    stop("Maximum likelihood failed to converge at initial step")
  scale <- fit$par.ests[1]
  shape <- fit$par.ests[2]
  theta <- c(scale, shape)
  if(bootstrap == FALSE & shape > 1)
    stop("Estimated parameters are outside the table range, please use the bootstrap version")
  thresh <- findthresh(data, n)
  newdata <- pgpd(data, loc = thresh, scale = scale, shape = shape)
  newdata <- sort(newdata)
  i <- seq(1, n, 1)
  stat <- -n - (1/n)*sum((2*i - 1)*(log(newdata) + log1p(-rev(newdata))))
  if(bootstrap == TRUE) {
    if(allowParallel == TRUE) {
      cl <- makeCluster(numCores)
      fun <- function(cl) {
        parSapply(cl, 1:bootnum, function(i,...) {gpdAdGen(n, theta)})
      }
      teststat <- fun(cl)
      stopCluster(cl)
    } else {
      teststat <- replicate(bootnum, gpdAdGen(n, theta))
    }
    teststat <- teststat[!is.na(teststat)]
    eff <- length(teststat)
    p <- (sum(teststat > stat) + 1) / (eff + 2)
  } else {
    row <- which(rownames(ADQuantiles) == max(round(shape, 2), -0.5))
    if(stat > ADQuantiles[row, 999]) {
      pvals <- -log(as.numeric(colnames(ADQuantiles[950:999])))
      x <- as.numeric(ADQuantiles[row, 950:999])
      y <- lm(pvals ~ x)
      stat <- as.data.frame(stat)
      colnames(stat) <- c("x")
      p <- as.numeric(exp(-predict(y, stat)))
    } else {
      bound <- as.numeric(colnames(ADQuantiles)[which.max(stat < ADQuantiles[row,])])
      if(bound == .999) {
        p <- .999
      } else {
        lower <- ADQuantiles[row, which(colnames(ADQuantiles) == bound + 0.001)]
        upper <- ADQuantiles[row, which(colnames(ADQuantiles) == bound)]
        dif <- (upper - stat) / (upper - lower)
        val <- (dif * (-log(bound) - -log(bound + 0.001))) + log(bound)
        p <- exp(val)
      }
    }
  }
  names(theta) <- c("Scale", "Shape")
  if(!bootstrap) {
    out <- list(as.numeric(stat), as.numeric(p), theta)
    names(out) <- c("statistic", "p.value", "theta")
  } else {
    out <- list(as.numeric(stat), as.numeric(p), theta, eff)
    names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  }
  out
}
