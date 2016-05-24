#' For given scores \code{Y} and dynamic principal components \code{XI}
#' retrive a series from which scores \code{Y} were calculated.
#' This procedure should be seen as the inverse of \code{\link{dpca.scores}}.
#'
#' @title Retrieve a process from given scores
#' @param Y scores process
#' @param XI principal components series
#' @return Retrived process X
#' @seealso \code{\link{dpca.scores}}, \code{\link{dprcomp}}
#' @references Siegfried Hormann, Lukasz Kidzinski and Marc Hallin
#' Dynamic Functional Principal Component
#' Research report, 2012
# @export
dpca.inverse = function(Y,XI){
  t(rev(XI)) %c% Y
}
