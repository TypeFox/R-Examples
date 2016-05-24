#' Blow up single scalars / objects to vectors /  list by replication.
#'
#' Useful for standard argument conversion where a user can input a single
#' element, but this has to be replicated now n times for a resulting vector or list.
#'
#' @param x [any]\cr
#'   Input element.
#' @param n [\code{integer}]\cr
#'   Desired length.
#' @param cl [\code{character(1)}*]\cr
#'   Only do the operation if \code{x} inherits from this class, otherwise simply let x pass.
#'   Default is \code{NULL} which means to always do the operation.
#' @param names [\code{character}*] \cr
#'   Names for result.
#'   Default is \code{NULL}, which means no names.
#' @return Ether a vector or list of length \code{n} with replicated \code{x} or \code{x} unchanged..
#' @export
ensureVector = function(x, n, cl = NULL, names = NULL) {
  n = convertInteger(n)
  assertCount(n)

  doit = isScalarValue(x) || !is.atomic(x)
  if (!is.null(cl)) {
    assertString(cl)
    doit = doit && inherits(x, cl)
  }

  if (doit) {
    if (isScalarValue(x))
      xs = rep(x, n)
    else
      xs = replicate(n, x, simplify = FALSE)

    if (!is.null(names)) {
      assertCharacter(names, len = n, any.missing = FALSE)
      names(xs) = names
    }
    return(xs)
  } else {
    return(x)
  }
}

