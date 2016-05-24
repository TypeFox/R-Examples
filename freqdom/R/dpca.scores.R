#' Compute scores of a given series \code{X} using dynamic principal components \code{XI}.
#' Procedure can be inverted using \code{\link{dpca.inverse}}
#'
#' @title Compute scores of dynamic principal components
#' @param X series to 'project'
#' @param XI principal components series
#' @return Matrix of scores of \code{X}
#' @seealso \code{\link{dpca.inverse}}, \code{\link{dprcomp}}
#' @references Siegfried Hormann, Lukasz Kidzinski and Marc Hallin
#' Dynamic Functional Principal Component
#' Research report, 2012
# @export
dpca.scores = function(X,XI){
  XI %c% X
}

