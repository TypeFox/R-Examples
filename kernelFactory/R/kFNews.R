#' Display the NEWS file
#'
#' \code{kFNews} shows the NEWS file of the kernelFactory package.
#' 
#' @examples
#' kFNews()
#' @return None.
#' @references Ballings, M. and Van den Poel, D. (2013), Kernel Factory: An Ensemble of Kernel Machines. Expert Systems With Applications, 40(8), 2904-2913.
#' @seealso \code{\link{kernelFactory}}, \code{\link{predict.kernelFactory}}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
kFNews <-
function() {
    newsfile <- file.path(system.file(package="kernelFactory"), "NEWS")
    file.show(newsfile)
}
