#' Pretty Breakpoints on Log Scale
#'
#' Compute a sequence of "round" values which cover the range of \code{x} on
#' the log scale.
#' @param x
#'   A numeric vector.
#' @param lead
#'   An integer vector giving the desired lead digits of pretty values on
#'   the log scale, default c(1, 5).
#' @param extra
#'   An integer scalar giving the desired number of additional
#'   non-log scale values to include, default 5.
#' @return
#'   A numeric vector of pretty values covering the range of \code{x} on
#'   the log scale.
#' @export
#' @examples
#' vals <- rlnorm(100, 6)
#' summary(vals)
#' prettylog(vals, 1, 0)
#' prettylog(vals, 1)
#' prettylog(vals, c(1, 2, 5))

prettylog <- function(x, lead=c(1, 5), extra=5) {
  if (!is.numeric(x)) stop("x must be a numeric vector")
  if (is.character(lead) | sum(abs(as.integer(lead) - lead))>0) {
    stop("lead must be an integer vector")
  }
  if (is.character(extra) | abs(as.integer(extra) - extra)>0 |
      length(extra)!=1) stop("extra must be an integer scalar")
  urd <- function(d, x) {
    lxd <- log10(x/d)
    rlxd <- unique(c(floor(lxd), ceiling(lxd)))
    d*10^rlxd
  }
  out <- sort(unlist(lapply(lead, urd, x)))
  if (extra>0) {
    out <- sort(unique(c(out, pretty(x, n=extra))))
  }
  out[out>0]
}
