#' Get Names Attributes of Growth Models
#'
#' Methods to get the parameter names of a growth model or to get or set
#'   identifiers of \code{\link{multiple_fits}} objects.
#'
#' @param x either a function being a parametric growth model of
#'   package \pkg{growthmodels} or an object with multiple fits.
#' @param value a character vector of up to the same length as x, or NULL
#'
#' @return character vector of the parameter names
#'
#' @section Methods:
#' \describe{
#'   \item{Method for class \code{\link{growthmodel}}:}{ returns information about
#'   valid parameter names if a \code{pnames} attribute exists, else \code{NULL}.
#'   \code{NULL}.}
#'   \item{Method for class \code{\link{multiple_fits}}:}{ can be applied to objects
#'   returned by \code{all_growthmodels}, \code{all_splines} or
#'   \code{all_easylinear} respectively. This can be useful for selecting
#'   subsets, e.g. for plotting, see example below.}
#' }

#'
#' @seealso \code{\link{multiple_fits}}, \code{\link{all_growthmodels}},
#'   \code{\link{all_splines}}, \code{\link{all_easylinear}}
#'
#' @examples
#'
#' ## growthmodel-method
#' names(grow_baranyi)
#'
#' ## multiple_fits-method
#' L <- all_splines(value ~ time | strain + conc + replicate,
#'        data = bactgrowth)
#'
#' names(L)
#'
#' ## plot only the 'R' strain
#' par(mfrow=c(4, 6))
#' plot(L[grep("R:", names(L))])
#'
#'
#' @rdname names
#' @export
#'
names.growthmodel <- function(x) attr(x, "pnames")

## S4 method does not work here (even if a setGeneric 'pnames' would do)
##    so we use the S3 method above
#setMethod("names", "growthmodel",
#          function(x) {
#            attr(x, "pnames")
#          }
#)


#' @rdname names
#' @exportMethod names
#'
setMethod("names", "multiple_fits",
          function(x) {
            names(x@fits)
          }
)

#' @rdname names
#' @exportMethod names<-
#'
setMethod("names<-", c("multiple_fits", "ANY"),
          function(x, value) {
            if (!is.character(value))
              value <- as.character(value)
            ## todo: check length?
            names(x@fits) <- value
          }
)


