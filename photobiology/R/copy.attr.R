#' Copy attributes
#'
#' Copy attributes from one R object to another.
#'
#' @param x,y R objects
#' @param which character
#' @param ... not used
#'
#' @return A copy of \code{y} with additional attributes set.
#'
#' @export
#'
copy_attributes <- function(x, y, which, ...) UseMethod("copy_attributes")

#' @describeIn copy_attributes Default for generic function
#'
#' @export
#'
copy_attributes.default <- function(x, y,
                                    which = NULL,
                                    ...) {
  warning("'copy_attributes' is not defined for objects of class ", class(x)[1])
  y
}

#' @describeIn copy_attributes
#'
#' @export
#'
copy_attributes.generic_spct <- function(x, y,
                                         which = NULL,
                                         ...) {
  stopifnot(is.any_spct(y))
  if (length(which) == 0) {
    which <- c("comment",
               "when.measured",
               "where.measured",
               "what.measured",
               "spct.tags",
               "normalized",
               "scaled",
               "multiple.wl",
               "spct.version")
    which.add <- c(switch, class(y),
                   generic_spct,
                   raw_spct,
                   cps_spct = c("instr.desc", "instr.settings"),
                   source_spct,
                   response_spct = c("time.unit", "bswf.used"),
                   # need to be copied in case class of object_spct
                   # is changed temporarily
                   filter_spct,
                   reflector_spct,
                   object_spct = c("Tfr.type", "Rfr.type"),
                   chroma_spct = character()
    )
    which <- c(which, which.add)
  }
  attr.x <- attributes(x)
  which.x <- intersect(names(attr.x), which)
  # this is likely to be slow
  for (w in which.x) {
        attr(y, w) <- attr(x, w)
  }
  y
}

#' @describeIn copy_attributes
#'
#' @export
#'
copy_attributes.waveband <- function(x, y, which = NULL, ...) {
  stopifnot(is.waveband(y))
  if (length(which) == 0) {
    which <- "comment"
  }
  attr.x <- attributes(x)
  which.x <- intersect(names(attr.x), which)
  for (w in which.x) {
    attr(y, w) <- attr(x, w)
  }
  y
}

