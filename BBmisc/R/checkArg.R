#' Check for a function argument.
#'
#' Throws exception if checks are not passed.
#' Note that argument is evaluated when checked.
#'
#' @param x [any]\cr
#'   Argument.
#' @param cl [\code{character}]\cr
#'   Class that argument must \dQuote{inherit} from.
#'   If multiple classes are given, \code{x} must \dQuote{inherit} from at least one of these.
#'   See also argument \code{s4}.
#' @param s4 [\code{logical(1)}]\cr
#'   If \code{TRUE}, use \code{is} for checking class \code{cl}, otherwise use \code{\link{inherits}}, which
#'   implies that only S3 classes are correctly checked. This is done for speed reasons
#'   as calling \code{\link{is}} is pretty slow.
#'   Default is \code{FALSE}.
#' @param len [\code{integer(1)}]\cr
#'   Length that argument must have.
#'   Not checked if not passed, which is the default.
#' @param min.len [\code{integer(1)}]\cr
#'   Minimal length that argument must have.
#'   Not checked if not passed, which is the default.
#' @param max.len [\code{integer(1)}]\cr
#'   Maximal length that argument must have.
#'   Not checked if not passed, which is the default.
#' @param choices [any]\cr
#'   Discrete number of choices, expressed by a vector of R objects.
#'   If passed, argument must be identical to one of these and nothing else is checked.
#' @param subset [any]\cr
#'   Discrete number of choices, expressed by a vector of R objects.
#'   If passed, argument must be identical to a subset of these and nothing else is checked.
#' @param lower [\code{numeric(1)}]\cr
#'   Lower bound for numeric vector arguments.
#'   Default is \code{NA}, which means not required.
#' @param upper [\code{numeric(1)}]\cr
#'   Upper bound for numeric vector arguments.
#'   Default is \code{NA}, which means not required.
#' @param na.ok [\code{logical(1)}]\cr
#'   Is it ok if a vector argument contains NAs?
#'   Default is \code{TRUE}.
#' @param formals [\code{character}]\cr
#'   If this is passed, \code{x} must be a function.
#'   It is then checked that \code{formals} are the names of the
#'   (first) formal arguments in the signature of \code{x}.
#'   Meaning \code{checkArg(function(a, b), formals = "a")} is ok.
#'   Default is missing.
#' @return Nothing.
#' @export
#' @examples
#' x = 1L
#' checkArg(x, "integer", len = 1, na.ok = FALSE, upper = 3L)
#' x = as.integer(NA)
#' checkArg(x, "integer", len = 1, na.ok = TRUE)
#' x = c("foo", "bar")
#' checkArg(x, "character")
#' x = "foo"
#' checkArg(x, choices = c("foo", "bar"))
#' x = c("foo", "bar")
#' checkArg(x, subset = c("foo", "bar"))
#' fun = function(foo, bar)
#' checkArg(fun, formals = c("foo", "bar"))
checkArg = function(x, cl, s4 = FALSE, len, min.len, max.len, choices, subset, lower = NA, upper = NA, na.ok = TRUE, formals) {
  s = deparse(substitute(x))
  if (missing(x))
    stop("Argument ", s, " must not be missing!")
  cl2 = class(x)[1]
  len2 = length(x)
  matchEl = function(x, xs) any(sapply(xs, function(y) identical(y, x)))
  # choices must be done first
  if (!missing(choices)) {
    if (!matchEl(x, choices))
      stop("Argument ", s, " must be any of: ", collapse(choices), "!")
  } else if (!missing(subset)) {
    if (!all(sapply(x, matchEl, xs = subset)))
      stop("Argument ", s, " must be subset of: ", collapse(subset), "!")
  } else if (!missing(formals)) {
    if (!is.function(x))
      stop("Argument ", s, " must be of class ", "function", " not: ", cl2, "!")
    fs = names(formals(x))
    if (length(fs) < length(formals) || !all(formals == fs[seq_along(formals)]))
      stop("Argument function must have first formal args: ", paste(formals, collapse = ","), "!")
  } else {
    mycheck = function(x, cc)
      if(identical(cc, "numeric"))
        is.numeric(x)
      else if(identical(cc, "integer"))
        is.integer(x)
      else if(identical(cc, "vector"))
        is.vector(x)
      else if (!s4)
        inherits(x, cc)
      else if (s4)
        is(x, cc)
    if (!any(sapply(cl, mycheck, x = x)))
      stop("Argument ", s, " must be of class ", collapse(cl, " OR "), ", not: ", cl2, "!")
    if (!missing(len) && len2 != len)
      stop("Argument ", s, " must be of length ", len, " not: ", len2, "!")
    if (!missing(min.len) && len2 < min.len)
      stop("Argument ", s, " must be at least of length ", min.len, " not: ", len2, "!")
    if (!missing(max.len) && len2 > max.len)
      stop("Argument ", s, " must be at most of length ", max.len, " not: ", len2, "!")
    if (!na.ok && any(is.na(x)))
      stop("Argument ", s, " must not contain any NAs!")
    if (is.numeric(x) && !is.na(lower) && ((is.na(x) && !na.ok) || (!is.na(x) && any(x < lower))))
      stop("Argument ", s, " must be greater than or equal ", lower, "!")
    if (is.numeric(x) && !is.na(upper) && ((is.na(x) && !na.ok) || (!is.na(x) && any(x > upper))))
      stop("Argument ", s, " must be less than or equal ", upper, "!")
  }
}
