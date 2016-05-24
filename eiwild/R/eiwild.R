#' eiwild
#'
#' @description
#' eiwild means \emph{\strong{E}cological \strong{I}nference \strong{W}ith \strong{I}ndividual 
#' \strong{L}evel \strong{D}ata}.
#' 
#' It enables to estimate the inner cells of an RxC-table with aggregate data with 
#' a number of precincts and/or individual level data for a few precincts.
#' 
#' The assumptions are based on \emph{Rosen et al's} (2001) on the Multinomial-Dirichlet-
#' Model and the \emph{Wakefield's} (2004) paper on combining individual level data for
#' 2x2-tables. In the master thesis of Thomas Schlesinger (2013) he expanded the 2x2-case
#' to the RxC-case and implemented it in this package. Some of the functions are 
#' based on the equivalent functions of the \emph{eiPack}-package.
#' 
#' It is a hierarchical Bayesian model which uses MCMC-Algorithms to calculate the
#' estimations. Therefore it has all the pros and cons of these methods.
#' 
#' The typical workflow consists of 
#' \enumerate{
#'   \item Formatting your data in the accepted way. 
#'      Loading \code{\link[eiwild]{topleveldat}} and inspecting the \code{data.frame}'s 
#'      \code{aggr} and \code{indi} will give you a general idea.
#'   \item Running one or multiple times:
#'   \enumerate{
#'     \item function \code{\link[eiwild]{tuneVars}}
#'     \item function \code{\link[eiwild]{indAggEi}}
#'   }   
#'   \item Analyzing results with \code{\link[eiwild]{plot.eiwild}},
#'      \code{\link[eiwild]{comPlot}}, \code{\link[eiwild]{summary.eiwild}},
#'      \code{\link[eiwild]{plotResult}}
#' }
#' 
#' For more examples look at the appropriate functions
#' 
#' @references
#' <<to be added>>
#' @name eiwild
#' @docType package
#' @aliases
#' eiwild
#' package-eiwild
#' 
#' @useDynLib eiwild
#' @import gtools
#' @import coda
#' @import lattice
#' 
#' @author Thomas Schlesinger 
NULL