#' @title Returns unique state assignment from a (row-wise) weight matrix 
#'
#' @description 
#' Converts a probabilistic cluster assignment to a unique cluster assignment
#' using the 
#' \describe{
#'   \item{\code{'argmax'} rule:}{state of row \eqn{i} is assigned as the
#'                                position of the maximum in that row 
#'                                (ties are broken at random).}
#'   \item{\code{'sample'} rule}{state of row \eqn{i} is sampled from the
#'                               discrete distribution where probabilities equal
#'                               the weight vector in row \eqn{i}}
#' }
#' 
#' 
#' @param weight.matrix an \eqn{N \times K} matrix
#' @param rule how do we choose the state given the weight matrix. 
#' \code{c("argmax", "sample")}.
#' @keywords manip array
#' @export
#' @seealso \code{\link{states2weight_matrix}}
#' @examples
#' WW = matrix(runif(12), ncol = 3)
#' WW = normalize(WW)
#' WW
#' weight_matrix2states(WW)
#' weight_matrix2states(WW, "sample")
#' # another 'sample' is in general different from previous conversion unless
#' # WW is a 0/1 matrix
#' weight_matrix2states(WW, "sample") 

weight_matrix2states <- function(weight.matrix, 
                                 rule = c("argmax", "sample")) {
  rule <- match.arg(rule)
  switch(rule,
         argmax = {
           states <- max.col(weight.matrix) # much faster than 'apply(weight.matrix, 1, which.max)'
         },
         sample = {
           # TODO: make this faster?
           states <- apply(weight.matrix, 1, 
                           function(pp) {
                             sample.int(ncol(weight.matrix), size = 1, prob = pp, 
                                        replace = TRUE)
                             })
         })
  invisible(states)
}
