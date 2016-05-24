#' Display the NEWS file
#'
#' \code{dummyNews} shows the NEWS file of the \code{dummy} package.
#' 
#' @examples
#' dummyNews()
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
dummyNews <-
function() {
    newsfile <- file.path(system.file(package="dummy"), "NEWS")
    file.show(newsfile)
}
