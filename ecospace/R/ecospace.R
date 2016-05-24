#' ecospace: Simulating Community Assembly and Ecological Diversification Using
#' Ecospace Frameworks
#'
#' ecospace is an R package that implements stochastic simulations of community
#' assembly (ecological diversification) using customizable ecospace frameworks
#' (functional trait spaces). Simulations model the 'neutral', 'redundancy',
#' 'partitioning", and "expansion" models of Bush and Novack-Gottshall (2012). It
#' provides a wrapper to calculate common ecological disparity and functional
#' ecology statistical dynamics as a function of species richness. Functions are
#' written so they will work in a parallel-computing environment.
#'
#' The package also contains a sample data set, functional traits for Late
#' Ordovician (Type Cincinnatian) fossil species from the Kope and Waynesville
#' formations, from Novack-Gottshall (In pressB).
#'
#' @author Phil Novack-Gottshall \email{pnovack-gottshall@@ben.edu}
#' @name ecospace-package
#' @aliases ecospace-package ecospace
#' @docType package
#'
#' @references Bush, A. and P.M. Novack-Gottshall. 2012. Modelling the
#' ecological-functional diversification of marine Metazoa on geological time
#' scales. \emph{Biology Letters} 8: 151-155.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology}, submitted
#' Oct. 5, 2015. General models of ecological diversification. I. Conceptual
#' synthesis.
#' @references Novack-Gottshall, P.M. In review at \emph{Paleobiology}, submitted
#' Oct. 5, 2015. General models of ecological diversification. II. Simulations
#' and empirical applications.
#'
#' @seealso The 'calc_metrics' function relies extensively on the functional
#' diversity package \code{\link[FD]{FD}}, and hence lists this package as a
#' depends, so it is loaded simultaneously.
#'
#' @examples
#' # Get the package version and citation of ecospace
#' packageVersion("ecospace")
#' citation("ecospace")
#'
#' # Create an ecospace framework (functional-trait space) with 15 characters
#' #    (functional traits) of mixed types
#' nchar <- 15
#' ecospace <- create_ecospace(nchar=nchar, char.state=rep(3, nchar),
#'   char.type=rep(c("factor", "ord.fac", "ord.num"), nchar / 3))
#'
#' # Use to assemble a stochastic "neutral" sample of 20 species (from
#' #    initial seeding by 5 species)
#' x <- neutral(Sseed=5, Smax=20, ecospace=ecospace)
#' head(x, 10)
#'
#' # Calculate ecological disparity (functional diversity) dynamics as a
#' #    function of species richness
#' # Statistic 'V' [total variance] not calculated because there are factors
#' #    in the sample
#' metrics <- calc_metrics(samples=x, Smax=20, Model="Neutral", Param="NA")
#' metrics
#'
#' # Plot statistical dynamics as function of species richness
#' op <- par()
#' par(mfrow=c(2,4), mar=c(4, 4, 1, .3))
#' attach(metrics)
#' plot(S, H, type="l", cex=.5)
#' plot(S, D, type="l", cex=.5)
#' plot(S, M, type="l", cex=.5)
#' plot(S, V, type="l", cex=.5)
#' plot(S, FRic, type="l", cex=.5)
#' plot(S, FEve, type="l", cex=.5)
#' plot(S, FDiv, type="l", cex=.5)
#' plot(S, FDis, type="l", cex=.5)
#'
#' par(op)
#'
# NAMESPACE IMPORTING
#' @import FD
#' @importFrom stats var rmultinom runif
NULL
