#' Display the NEWS file
#'
#' \code{hybridEnsembleNews} shows the NEWS file of the hybridEnsemble package.
#' 
#' @examples
#' hybridEnsembleNews()
#' @return None.
#' @references Ballings, M., Vercamer, D., Van den Poel, D., Hybrid Ensemble: Many Ensembles is Better Than One, Forthcoming.
#' @seealso \code{\link{hybridEnsemble}}, \code{\link{predict.hybridEnsemble}}, \code{\link{importance.hybridEnsemble}}, \code{\link{CVhybridEnsemble}}, \code{\link{plot.CVhybridEnsemble}}, \code{\link{summary.CVhybridEnsemble}}
#' @author Michel Ballings, Dauwe Vercamer, and Dirk Van den Poel, Maintainer: \email{Michel.Ballings@@GMail.com}
hybridEnsembleNews <-
function() {
    newsfile <- file.path(system.file(package="hybridEnsemble"), "NEWS")
    file.show(newsfile)
}
