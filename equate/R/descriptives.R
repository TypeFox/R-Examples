#' Descriptive Statistics for Frequency Tables
#' 
#' These functions return descriptive statistics for a frequency table of class
#' \dQuote{\code{\link{freqtab}}}.
#' 
#' \code{mean}, \code{sd.freqtab}, \code{var.freqtab}, \code{skew.freqtab}, and
#' \code{kurt.freqtab} return the mean, standard deviation, variance, skewness,
#' and kurtosis. \code{min} and \code{max} return the minimum and maximum
#' observed scores, and \code{range} returns both. \code{cov.freqtab} and
#' \code{cor.freqtab} return the covariance and correlation matrices for one or
#' more variables. \code{summary} returns univariate statistics across one or
#' more margins.
#' 
#' @param object,x object of class \dQuote{\code{freqtab}}.
#' @param margin integer vector specifying the margin(s) for which summary
#' statistics will be returned. This defaults to \code{1} for univariate
#' statistics, and \code{seq(margins(x))}, i.e., all the margins, for
#' multivariate statistics (covariance and correlation).
#' @param \dots further arguments passed to or from other methods.
#' @param na.rm logical indicating whether missing values should be removed,
#' currently ignored since frequency tables cannot contain missing values.
#' @return \code{summary} returns a data frame of summary statistics, including
#' the mean, standard deviation, skewness, kurtosis, minimum, maximum, and
#' number of observations for each variable in \code{margin}. Otherwise, a
#' vector of length \code{length(margin)} is returned with the corresponding
#' statistic for each variable.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{freqtab}}
#' @keywords methods
#' @examples
#' 
#' summary(as.freqtab(ACTmath[, 1:2]))
#' 
#' ny <- freqtab(KBneat$y, scales = list(0:36, 0:12))
#' summary(ny)
#' cov.freqtab(ny)
#' @export
summary.freqtab <- function(object,
  margin = seq(margins(object)), ...) {
  
  out <- NULL
  for (i in margin) {
    xm <- margin(object, i)
    out <- rbind(out, data.frame(
      mean = mean(xm),
      sd = sd.freqtab(xm),
      skew = skew.freqtab(xm),
      kurt = kurt.freqtab(xm),
      min = min(xm),
      max = max(xm),
      n = sum(xm)))
  }
  rownames(out) <- names(dimnames(object))[margin]
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
mean.freqtab <- function(x, margin = 1, ...) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  out <- sapply(margin, function(i) {
    xm <- margin(x, i)
    sum(xm * scales(xm)/sum(xm))})
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
sd.freqtab <- function(x, margin = 1) {
  
  return(sqrt(var.freqtab(x, margin)))
}

#' @rdname summary.freqtab
#' @export
var.freqtab <- function(x, margin = 1) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  n <- sum(x)
  out <- sapply(margin, function(i) {
    xm <- margin(x, i)
    xsc <- scales(xm)
    (sum(xsc * xsc * xm) - (sum(xm * xsc)^2)/n)/(n - 1)})
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
cov.freqtab <- function(x, margin = seq(margins(x))) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  n <- sum(x)
  nx <- length(margin)
  out <- matrix(nrow = nx, ncol = nx)
  for (i in 1:nx) {
    out[i, i:nx] <- out[i:nx, i] <-
      sapply(margin[i:nx], function(j) {
        xd <- as.data.frame(margin(x,
          unique(c(i, j))))
        nc <- ncol(xd)
        sum((xd[, 1] - mean(x, i)) *
            (xd[, nc - 1] - mean(x, j)) *
            xd[, nc])/(n - 1)})
  }
  
  attr(out, "dim") <- c(nx, nx)
  attr(out, "dimnames") <- list(names(dimnames(x))[margin],
    names(dimnames(x))[margin])
  return(out)
}

#' @rdname summary.freqtab
#' @export
cor.freqtab <- function(x, margin = seq(margins(x))) {
  
  sds <- 1/sd.freqtab(x, margin)
  covs <- cov.freqtab(x, margin)
  out <- diag(sds) %*% covs %*% diag(sds)
  attributes(out) <- attributes(covs)
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
min.freqtab <- function(x, margin = 1, ..., na.rm = FALSE) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  out <- sapply(margin, function(i) {
    xm <- margin(x, i)
    min(scales(xm)[as.logical(xm)])})
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
max.freqtab <- function(x, margin = 1, ..., na.rm = FALSE) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  out <- sapply(margin, function(i) {
    xm <- margin(x, i)
    max(scales(xm)[as.logical(xm)])})
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
range.freqtab <- function(x, margin = 1, ...,
  na.rm = FALSE) {
  
  if (length(margin) == 1)
    out <- c(min(x), max(x))
  else if (length(margin) > 1)
    out <- lapply(margin, function(i)
      c(min(margin(x, i)), max(margin(x, i))))
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
skew.freqtab <- function(x, margin = 1) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  n <- sum(x)
  out <- sapply(margin, function(i) {
    xm <- margin(x, i)
    xsc <- scales(xm)
    sum(((xsc - mean(xm))^3 * xm))/(n)/
      (sum(((xsc - mean(xm))^2 * xm))/(n - 1))^1.5})
  
  return(out)
}

#' @rdname summary.freqtab
#' @export
kurt.freqtab <- function(x, margin = 1) {
  
  inmars <- margin %in% seq(margins(x))
  if (any(!inmars)) {
    margin <- margin[inmars]
    warning("misspecified margins ",
      "have been removed.")
  }
  n <- sum(x)
  out <- sapply(margin, function(i) {
    xm <- margin(x, i)
    xsc <- scales(xm)
    sum(((xsc - mean(xm))^4 * xm))/(n)/
      (sum(((xsc - mean(xm))^2 * xm))/(n - 1))^2})
  
  return(out)
}
