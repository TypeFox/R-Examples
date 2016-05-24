
#' Coercion from other classes to \code{trip} objects
#'
#' Coercing objects to \code{trip} class
#'
#' @name as.trip
#' @aliases as.trip-methods as.trip as.trip,ltraj-method ltraj2trip
#' coerce,trip,ltraj-method
#' @docType methods
#' @param x, ltr ltraj object
#' @param \dots Arguments passed to other methods. Ignored for \code{ltraj}
#' method.
#' @section Methods:
#'
#' \describe{
#'
#' \item{coerce}{\code{signature(from="ltraj", to="trip")}}
#'
#' \item{as.trip}{\code{signature(x="ltraj")}}
#'
#' }
#' @examples
#' ## Continuing the example from '?trip-methods:
#' utils::example("trip-methods", package="trip",
#'                ask=FALSE, echo=FALSE)
#'
#' if (require(adehabitatLT)) {
#'     ##l <- as.ltraj.trip(tr)
#'     ##ltraj2trip(l)
#'     ##as.trip(l)
#' }
##' @rdname as.trip-methods
##' @export
setGeneric("as.trip",
           function(x, ...) standardGeneric("as.trip"))

##' @export
ltraj2trip <- function (ltr)
{
  requireNamespace("adehabitatLT") ||
    stop("adehabitatLT package is required, but unavailable")
  if (!inherits(ltr, "ltraj"))
    stop("ltr should be of class \"ltraj\"")
  ltr <-  lapply(ltr, function(x) {
    x$id=attr(x,  "id")
    x$burst=attr(x,  "burst")
    x})
  tr <- do.call("rbind", ltr)
  class(tr) <- "data.frame"
  xy <- tr[!is.na(tr$x), c("x", "y")]
  tr <- tr[!is.na(tr$x), ]
  tr$y <- tr$x <- NULL
  res <- SpatialPointsDataFrame(xy, tr)
  trip(res, c("date", "id"))
}

setMethod("as.trip", signature(x="ltraj"),
          function(x, ...) ltraj2trip(x))

setAs("ltraj", "trip", function(from) as.trip(from))

