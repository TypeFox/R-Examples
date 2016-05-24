#' Extraction and assignment of parameters to vector
#'
#' Assigns parameters of a particular type (scale,
#' shape, adjustments or g0 (p(0))) from the vector of parameters in
#' \code{ddfobj}. All of the parameters are kept in a single vector for
#' optimization even though they have very different uses.  \code{assign.par}
#' parses them from the vector based on a known structure and assigns them into
#' \code{ddfobj}.  \code{getpar} extracts the requested types to be extracted
#' from \code{ddfobj}.
#'
#' @aliases assign.par
#' @param ddfobj distance sampling object (see \code{\link{create.ddfobj}})
#' @param fpar parameter vector
#' @return index==FALSE, vector of parameters that were requested or
#'   index==TRUE, vector of 3 indices for scale, shape, adjustment
#' @note Internal functions not intended to be called by user.
#' @seealso getpar
#' @author Jeff Laake
#' @keywords utility
assign.par <- function(ddfobj, fpar){
  index <- 1
  # Shape parameters
  if(!is.null(ddfobj$shape)){
    ddfobj$shape$parameters <- fpar[1:ncol(ddfobj$shape$dm)]
    index <- index+ncol(ddfobj$shape$dm)
  }

  # Scale parameters
  if(!is.null(ddfobj$scale)){
    ddfobj$scale$parameters <- fpar[index:(index+ncol(ddfobj$scale$dm)-1)]
    index <- index+ncol(ddfobj$scale$dm)
  }

  # Adjustment parameters
  if(!is.null(ddfobj$adjustment)){
    ddfobj$adjustment$parameters <- fpar[index:length(fpar)]
  }

  return(ddfobj)
}
