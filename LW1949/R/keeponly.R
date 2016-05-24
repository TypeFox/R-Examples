#' Eliminate Consecutive Extreme Values
#'
#' Generate the index for eliminating values beyond a given maximum number of
#'   consecutive extremes allowed.
#' @param x
#'   A numeric vector, with no missing values.
#' @param extremes
#'   A numeric vector of length two giving the boundary limits for \code{x},
#'     default c(0, 100).
#' @param nconsec
#'   An integer scalar, the maximum number of consecutive extreme values
#'     allowed, default 2.
#' @return
#'   A logical vector for selecting all elements of \code{x} without
#'     exceeding \code{nconsec} consecutive extreme values.
#' @export
#' @examples
#' vec <- c(0, 0, 0, 4, 4, 4, 100, 100, 100, 100)
#' vec[keeponly(vec)]
#' # the original vector need not be ordered
#' vec <- c(100, 4, 100, 4, 0, 100, 0, 4, 0, 100)
#' keeponly(vec)

keeponly <- function(x, extremes=c(0, 100), nconsec=2) {
  if (any(is.na(x))) stop("No missing values allowed in x.")
  ord <- order(x)
  orderedx <- x[ord]
  hi <- sapply(orderedx, all.equal, extremes[2])=="TRUE"
  selnohi <- cumsum(hi) <= nconsec
  lo <- rev(sapply(orderedx, all.equal, extremes[1])=="TRUE")
  selnolo <- rev(cumsum(lo) <= nconsec)
  ko <- selnolo & selnohi
  ko[order(ord)]
  }
