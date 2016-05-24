#' Display the NEWS file
#'
#' \code{interpretRNews} shows the NEWS file of the interpretR package.
#' 
#' @examples
#' interpretRNews()
#' @seealso \code{\link{parDepPlot}}
#' @author Authors: Michel Ballings and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
interpretRNews <-
function() {
    newsfile <- file.path(system.file(package="interpretR"), "NEWS")
    file.show(newsfile)
}
