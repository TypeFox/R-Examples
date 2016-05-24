#' @title Print the MEW average to the screen
#'
#' @description Print to the screen the first six elements of the
#' current value (if it is up-to-date) of the moving expanding window
#' (MEW) average. An error is raised if the MEW average is not
#' up-to-date.
#'
#' @return Upon successful exit, zero is returned invisibly.
#'
#' @param object (class mewTyp) The current state of the MEW average
#'
#' @examples
#' ## see the examples for the function mewMean
#'
#' @import methods
#'
#' @export
setMethod("show",
          signature = signature("mewTyp"),
          definition = function (object) {

            if(object@know_mean == 0) {

              stop("The mean is not up-to-date")
            } else {

              cat("Only the first six elements are shown\n")
              print(head(object@x_mean))
            }

            invisible(0)
          })
