#' @title Estimates LSC from a LICORS estimate
#'
#' @description 
#' A wrapper of \code{\link{states2LSC}} for a \code{'LICORS'} estimate from the
#' \pkg{LICORS} package (in particular the output from 
#' the \code{\link[LICORS]{mixed_LICORS}} 
#' function). Converts LICORS estimates into an array with LSC.
#' 
#' @param object an object of class \code{"LICORS"}
#' @param type should marginal state probabilities be computed based on the 
#' unique state space assignment (\code{"argmax"}) or using the soft 
#' thresholding from the weight matrix (\code{"weights"}).
#' @keywords manip
#' @return
#' An object of class \code{"LSC"}
#' @export
#' @seealso \code{\link{states2LSC}}, \code{\link{LSC-utils}}
#' @examples
#' \dontrun{
#' # see 2nd example in 'LSC-package'
#' }

LICORS2LSC = function(object, type = c("weights", "argmax")) {
  
  type <- match.arg(type)
  switch(type,
         argmax = {
           out <- states2LSC(states = object$states)
         },
         weights = {
           out <- states2LSC(weight.matrix = object$conditional.state.probs$opt)
         })

  space.dim <- length(object$dim$truncated) - 1
  
  if (space.dim == 0) {
    out <- ts(out)
  } else if (space.dim == 1) {
    out <- t(matrix(out, ncol = object$dim$truncated[2], byrow = TRUE))
  } else if (space.dim == 2) {
    out <- array(out, dim = c(object$dim$truncated[-1], object$dim$truncated[1]))
    # R.utils method
    #out = wrap( array(out, dim = c(object$dim$truncated[-1], object$dim$truncated[1])), map = list(2,1) )
  }

  class(out) <- "LSC"
  invisible(out)
}
