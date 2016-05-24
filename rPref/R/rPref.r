#' Summary of the rPref Package
#' 
#' rPref contains routines to select and visualize the maxima for a given strict
#' partial order. This especially includes the computation of the Pareto
#' frontier, also known as (Top-k) Skyline operator, and some
#' generalizations (database preferences).
#' 
#' @section Preference Composition/Selection:
#' 
#' \itemize{
#' \item Preferences are primarily composed from base preferences (see \code{\link{base_pref}}) 
#'   and complex preferences (see \code{\link{complex_pref}}),
#'   where especially the Pareto operator for Skylines is such a complex preference. 
#' \item Some utility functions for preferences are collected in \code{\link{general_pref}}.
#' \item Additionally some base preference macros are provided in \code{\link{base_pref_macros}}. 

#' \item The (top(-level)-k) preference selection \code{\link{psel}} allows to retrieve 
#'       the maxima of a preference (or Pareto frontier, Skyline), 
#'       constructed with the functions above, on a given data set.
#' }
#' 
#' @section Visualization and Analysis of Preferences:
#' 
#' \itemize{
#' \item The visualization of the preference order in a Better-Than-Graph (Hasse diagram) is possible via the \code{\link{plot_btg}} 
#'       and \code{\link{get_btg}} functions in connection with the \code{\link{igraph}} package.
#' \item The adjacency list of the Hasse diagram can be accessed via \code{\link{get_hasse_diag}}.
#' \item Predecessors/successors in the Hasse diagram are calculated with the \code{\link{pred_succ}} functions.
#' \item The Pareto frontier can be plotted using the \code{\link{plot_front}} function.
#' }
#' 
#' @section String Output of Preferences:
#' 
#' \itemize{
#' \item The preference query for some preference-supporting DBMS can be given by \code{\link{show.query}}.
#' \item A preference is partially evaluated and printed with \code{\link{show.pref}}.
#' }
#' 
#' @section Further Information:
#' 
#' The rPref homepage is \url{http://www.p-roocks.de/rpref}. To submit bugs, feature requests or other comments, feel free to write a mail to me.
#' 
#' @author Patrick Roocks, \email{mail@@p-roocks.de}
#'
#' @docType package
#' @name rPref
#' @useDynLib rPref
#' @importFrom Rcpp cppFunction
#' @importFrom RcppParallel RcppParallelLibs
#' @import igraph methods
NULL
