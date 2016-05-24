#' Simulating the evolution of biological sequences
#'
#' The Sequential Coalescent with Recombination Model (SCRM) is an
#' approximation of the Ancestral Recombination Graph. It can be used to
#' simulate the neutral evolution of chromosomes/biological sequences
#' subject to possibly complicated population structure.
#' The program \emph{scrm} is an implementation of this model that is designed to
#' act as an drop-in replacement for the widely adopted coalescent simulator
#' \emph{ms}. This package contains \emph{scrm} along with a convinient R interface.
#'
#' @author
#' Paul Staab,
#' Zhu Sha,
#' Dirk Metzler &
#' Gerton Lunter
#'
#' Maintainer: Paul Staab \email{develop@@paulstaab.de}
#' @name scrm-package
#' @docType package
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib scrm
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{scrm}} for details on how to use \emph{scrm},
#'   \item \code{vignette('scrm-Arguments')} for an overview of command line arguments and
#'   \item \code{vignette('scrm-TreesForApe')} for an example on using
#'         genealogies simulated with \emph{scrm} with package 'ape'.
#' }
NULL
