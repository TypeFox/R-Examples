#' @title Checks if values are outside of an interval (open on the
#'   right)
#' @description This function checks if the values in the \code{x}
#'   parameter are contained in the interval [\code{lo}, \code{hi}).
#'   \code{NA} values are treated as "not in the interval".
#' @param x A vector of values. (Lists will be coerced to a numeric vector.)
#' @param lo Left end of the interval.
#' @param hi Right end of the interval.
#' @return A boolean vector of the same length as \code{x}.
#' @seealso \link{in.interval.lo}, \link{in.interval.ro},
#'   \link{nin.interval.lo}
#' @examples
#' nin.interval.ro(c(-1, 0, 1, 2), 0, 1)
#' nin.interval.ro(NA, 1, 3)
#' @export
nin.interval.ro <- function(x, lo, hi) is.na(.bincode(x, c(lo, hi), right=F))

#' @title Checks if values are contained in an interval (open on the
#'   right)
#' @description This function checks if the values in the \code{x}
#'   parameter are contained in the interval [\code{lo}, \code{hi}).
#'   \code{NA} values are treated as "not in the interval".
#' @param x A vector of values. (Lists will be coerced to a numeric
#'   vector.)
#' @param lo Left end of the interval.
#' @param hi Right end of the interval.
#' @return A boolean vector of the same length as \code{x}.
#' @seealso \link{in.interval.lo}, \link{nin.interval.lo},
#'   \link{nin.interval.ro}
#' @examples
#' in.interval.ro(c(-1, 0, 1, 2), 0, 1)
#' in.interval.ro(NA, 1, 3)
#' @export
in.interval.ro <- function(x, lo, hi) !nin.interval.ro(x, lo, hi)

#' @title Checks if values are outside of an interval (open on the left)
#' @description This function checks if the values in the \code{x}
#'   parameter are contained in the interval (\code{lo}, \code{hi}].
#'   \code{NA} values are treated as "not in the interval".
#' @param x A vector of values. (Lists will be coerced to a numeric
#'   vector.)
#' @param lo Left end of the interval.
#' @param hi Right end of the interval.
#' @return A boolean vector of the same length as \code{x}.
#' @seealso \link{in.interval.lo}, \link{in.interval.ro},
#'   \link{nin.interval.ro}
#' @examples
#' nin.interval.lo(c(-1, 0, 1, 2), 0, 1)
#' nin.interval.lo(NA, 1, 3)
#' @export
nin.interval.lo <- function(x, lo, hi) is.na(.bincode(x, c(lo, hi), right=T))

#' @title Checks if values are contained in an interval (open on the
#'   left)
#' @description This function checks if the values in the \code{x}
#'   parameter are contained in the interval (\code{lo}, \code{hi}].
#'   \code{NA} values are treated as "not in the interval".
#' @param x A vector of values. (Lists will be coerced to a numeric
#'   vector.)
#' @param lo Left end of the interval.
#' @param hi Right end of the interval.
#' @return A boolean vector of the same length as \code{x}.
#' @seealso \link{in.interval.ro}, \link{nin.interval.lo},
#'   \link{nin.interval.ro}
#' @examples
#' in.interval.lo(c(-1, 0, 1, 2), 0, 1)
#' in.interval.lo(NA, 1, 3)
#' @export
in.interval.lo <- function(x, lo, hi) !nin.interval.lo(x, lo, hi)
