# Get a vector of data depths for k letter values
getDepth <- function(k, n) {
  # compute letter values based on depth
  depth <- rep(0,k)
  depth[1] <- (1+n)/2

  if (k > 1) {
    for (j in 2:k) depth[j] <- (1 + floor(depth[j - 1])) / 2
  }
  depth
}

# Calculate first k letter values for vector x
calcLV <- function(x, k) {
  n <- length(x)
  depth <- getDepth(k, n)

  y <- sort(x)
  d <- c(rev(depth),n-depth+1)
  qu <- (y[floor(d)] + y[ceiling(d)])/2
  # floor and ceiling is the same for .0 values
  # .5 values yield average of two neighbours

  # k, k+1 is the median
  # report only once
  qu[-k]
}

# Determine names of first k letter values
# output is (1) list of names starting with 'M'edian, and
# (2) a vector of letter values ordered by rank from lower kth letter value to upper k letter value
nameLV <- function(k) {
  # list of
  #  k letter values (starting with median)
  #  lower/upper letter values ordered from lowest to highest
  lvs <- NULL
  conf <- "M"
  if (k > 1) {
    lvs <- c(LETTERS[6:1], LETTERS[c(26:14, 12:7)])[1:(k - 1)]
    conf <- c(paste(rev(lvs), "l", sep = ""), "M", paste(lvs, "u", sep = ""))
  }
  list(LV = c("M", lvs), conf = conf)
}

#' Determine depth of letter values needed for n observations.
#'
#' @details Supply one of \code{k}, \code{alpha} or \code{perc}.
#'
#' @param n number of observation to be shown in the LV boxplot
#' @param k number of letter value statistics used
#' @param alpha if supplied, depth k is calculated such that (1-\code{alpha})100% confidence
#'   intervals of an LV statistic do not extend into
#'   neighboring LV statistics.
#' @param perc if supplied, depth k is adjusted such that \code{perc} percent
#'   outliers are shown
#' @export
determineDepth <- function(n, k = NULL, alpha = NULL, perc = NULL) {
  if (!is.null(k)) {
    stopifnot(is.numeric(k) && length(k) == 1)
    k <- as.integer(k)
  } else if (!is.null(perc)) {
    # we're aiming for perc percent of outlying points
    stopifnot(is.numeric(perc) && length(perc) == 1)

    k <- ceiling((log2(n))+1) - ceiling((log2(n*perc*0.01))+1) + 1
  } else if (!is.null(alpha)) {
    # confidence intervals around an LV statistic
    # should not extend into surrounding LV statistics
    stopifnot(is.numeric(alpha) && length(alpha) == 1)
    stopifnot(alpha > 0 && alpha < 1)

#    cat(sprintf("two rules: %d (ceiling) %d (floor)", ceiling((log2(n))-log2(2*qnorm(alpha+(1-alpha)/2)^2)),
#            floor(log2(n)) - floor(log2(2*qnorm(1-(1-alpha)/2)^2))))
    k <- floor(log2(n)) - floor(log2(2*qnorm(1-(1-alpha)/2)^2))
  } else {
    stop("Must specify one of k, alpha, perc", call. = FALSE)
  }

  max(k, 1L)
}

#' Compute table of k letter values for vector x
#'
#' @param x input numeric vector
#' @param k number of letter values to compute
#' @param alpha alpha-threshold for confidence level
#' @export
lvtable <- function(x, k, alpha=0.95) {
  n <- length(x)
  if (2^k > n) k <- ceiling(log2(n)) + 1

  # depths for letter values
  depth <- getDepth(k, n)

  # letter value
  qu <- calcLV(x,k)

  tab <- matrix(c(c(rev(depth), depth[-1]), qu), ncol = 2,
    dimnames = list(nameLV(k)[[2]], c("depth","LV")))

  # confidence limits
  conf <- confintLV(x, k, alpha = alpha)

  cbind(tab, conf)
}


# confidence interval for k letter values
confintLV <- function(x, k, alpha=0.95) {
  n <- length(x)
  y <- sort(x)

  depth <- getDepth(k,n)
  extend <- ceiling(0.5 *sqrt(2*depth-1) * qnorm(alpha+(1-alpha)/2))
  low <- depth - extend
  high <- depth + extend
  clow <- pmax(1,ceiling(low))
  flow <- pmax(1,floor(low))
  chigh <- pmin(n, ceiling(high))
  fhigh <- pmin(n, floor(high))

  lvllow <- rev(rowMeans(cbind(y[clow],y[flow]), na.rm=T))
  if (length(lvllow) == 0) lvllow <- NA
  lvlhigh <- rev(rowMeans(cbind(y[chigh],y[fhigh]), na.rm=T))
  if (length(lvlhigh) == 0) lvlhigh <- NA
  # no 1 is the median - that's the last element in lvl
  lvulow <- rowMeans(cbind(y[n-chigh],y[n-fhigh]), na.rm=T)[-1]
  lvuhigh <- rowMeans(cbind(y[n-clow],y[n-flow]), na.rm=T)[-1]

  conf <- cbind(c(lvllow, lvulow), c(lvlhigh, lvuhigh))

  colnames(conf) <- c(paste((1 - alpha) / 2 * 100, "%" ,sep = "") ,
    paste((alpha + (1 - alpha) / 2) * 100, "%", sep = ""))
  conf
}
