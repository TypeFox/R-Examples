#' Iterators for sequence generation
#'
#' Constructs iterators that generate regular sequences that follow the
#' \code{\link[base]{seq}} family.
#'
#' The \code{iseq} function generates a sequence of values beginning with
#' \code{from} and ending with \code{to}. The sequence of values between are
#' determined by the \code{by}, \code{length_out}, and \code{along_with}
#' arguments. The \code{by} argument determines the step size of the sequence,
#' whereas \code{length_out} and \code{along_with} determine the length of the
#' sequence. If \code{by} is not given, then it is determined by either
#' \code{length_out} or \code{along_with}. By default, neither are given, in
#' which case \code{by} is set to 1 or -1, depending on whether \code{to >
#' from}.
#'
#' \code{seq_along} and \code{seq_len} return an iterator, which generates a
#' sequence of integers, beginning with 1 and proceeding to an ending value
#'
#' @export
#' @param from the starting value of the sequence
#' @param to the end value of the sequence
#' @param by increment of the sequence.
#' @param length_out desired length of the sequence. A non-negative number,
#' which for \code{seq} will be rounded up if fractional.
#' @param along_with the length of the sequence will match the length of this
#' argument
#' @return sequence's iterator
#' 
#' @examples
#' it <- iseq(from=2, to=5)
#' unlist(as.list(it)) == 2:5
#'
#' it2 <- iseq_len(4)
#' unlist(as.list(it2)) == 1:4
#'
#' it3 <- iseq_along(iris)
#' unlist(as.list(it3)) == 1:length(iris)
iseq <- function(from=1, to=1, by=(to - from)/(length_out - 1),
                 length_out=NULL, along_with=NULL) {
  from <- as.numeric(from)
  to <- as.numeric(to)
  by <- as.numeric(by)

  if (length(from) != 1) {
    stop("'from' must be a numeric value of length 1")
  }
  if (length(to) != 1) {
    stop("'to' must be a numeric value of length 1")
  }
  if (length(by) > 1 || (length(by) == 1 && by == 0)) {
    stop("'by' must be a nonzero numeric value of length 1")
  }

  # If 'by' is not given, then it is determined by either 'length_out' or
  # 'along_with'.  # By default, neither are given, in which case 'by' is set to
  # 1 or -1, depending on whether to > from.
  if (length(by) == 0) {
    if (!is.null(along_with)) {
      length_out <- length(along_with)
      by <- (to - from) / (length_out - 1)
    } else if (to >= from) {
      by <- 1
    } else if (to < from) {
      by <- -1
    }
  }

  current_val <- from - by
  nextElem <- function() {
    current_val <<- current_val + by
    if ((by > 0 && current_val > to) || ((by < 0) && current_val < to)) {
      stop("StopIteration", call.=FALSE)
    }
    current_val
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}

#' @export
#' @rdname iseq
iseq_len <- function(length_out=NULL) {
  length_out <- suppressWarnings(try(as.integer(length_out), silent=TRUE))

  if (inherits(length_out, "try-error") || length(length_out) != 1
      || length_out < 0 || is.na(length_out)) {
    stop("'length_out' must be coercible to non-negative integer")
  }

  i <- 0
  nextElem <- function() {
    i <<- i + 1

    if (i > length_out) {
      stop("StopIteration", call.=FALSE)
    }

    i
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}

#' @export
#' @rdname iseq
iseq_along <- function(along_with=NULL) {
  length_out <- length(along_with)

  i <- 0
  nextElem <- function() {
    i <<- i + 1

    if (i > length_out) {
      stop("StopIteration", call.=FALSE)
    }

    i
  }

  it <- list(nextElem=nextElem)
  class(it) <- c("abstractiter", "iter")
  it
}
