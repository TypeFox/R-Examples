#' Is given argument an atomic vector or factor of length 1?
#'
#' More specific functions for scalars of a given type exist, too.
#'
#' @param x [any]\cr
#'   Argument.
#' @param na.ok [\code{logical(1)}]\cr
#'   Is \code{NA} considered a scalar?
#'   Default is \code{TRUE}.
#' @param null.ok [\code{logical(1)}]\cr
#'   Is \code{NULL} considered a scalar?
#'   Default is \code{FALSE}.
#' @param type [\code{character(1)}]\cr
#'   Allows to restrict to specific type, e.g., \dQuote{numeric}?
#'   But instead of this argument you might want to consider using \code{isScalar<Type>}.
#'   Default is \dQuote{atomic}, so no special restriction.
#' @return [\code{logical(1)}].
#' @export
isScalarValue = function(x, na.ok = TRUE, null.ok = FALSE, type = "atomic") {
  if (is.null(x))
    return(null.ok)
  # not really cool switch, but maybe fastest option
  istype = switch(type,
    "atomic" = is.atomic,
    "logical" = is.logical,
    "numeric" = is.numeric,
    "integer" = is.integer,
    "complex" = is.complex,
    "chararacter" = is.character,
    "factor" = is.factor
  )
  istype(x) && length(x) == 1L && (na.ok || !is.na(x))
}

#' @rdname isScalarValue
#' @export
isScalarLogical = function(x, na.ok = TRUE, null.ok = FALSE) {
  isScalarValue(x, na.ok, null.ok, "logical")
}

#' @rdname isScalarValue
#' @export
isScalarNumeric = function(x, na.ok = TRUE, null.ok = FALSE) {
  isScalarValue(x, na.ok, null.ok, "numeric")
}

#' @rdname isScalarValue
#' @export
isScalarInteger = function(x, na.ok = TRUE, null.ok = FALSE) {
  isScalarValue(x, na.ok, null.ok, "integer")
}

#' @rdname isScalarValue
#' @export
isScalarComplex = function(x, na.ok = TRUE, null.ok = FALSE) {
  isScalarValue(x, na.ok, null.ok, "complex")
}

#' @rdname isScalarValue
#' @export
isScalarCharacter = function(x, na.ok = TRUE, null.ok = FALSE) {
  isScalarValue(x, na.ok, null.ok, "chararacter")
}

#' @rdname isScalarValue
#' @export
isScalarFactor = function(x, na.ok = TRUE, null.ok = FALSE) {
  isScalarValue(x, na.ok, null.ok, "factor")
}

