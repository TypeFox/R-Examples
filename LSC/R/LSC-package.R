#' @title
#' Local Statistical Compelxity - Automated Pattern Discovery in Spatio-Temporal Data
#' @aliases LSC
#' @docType package
#' @name LSC-package
#' @author
#' Georg M. Goerg \email{gmg@@stat.cmu.edu}
#' @description
#' A package to estimate local statistical complexity (LSC), a measure
#' for automated pattern discovery in spatio-temporal data using
#' optimal predictors (see References).
#' 
#' This package is very tightly linked to the \pkg{LICORS} package, which can
#' be used to estimate these optimal predictors and state space from data.  
#' The \pkg{LSC} builds on a known or estimated state space; most estimation
#' is handled by \pkg{LICORS} (see \code{?LICORS}).
#' 
#' There are two ways the state space can be represented: either as a unique
#' state label or as a vector of weights.  These two are the principal arguments 
#' in the functions of this package:
#' \describe{
#'  \item{\code{weight.matrix}}{an \eqn{N \times K} matrix, where \eqn{N} are 
#'                              the samples and \eqn{K} are the states. That is, 
#'                              each row contains a vector of length \eqn{K} that
#'                              adds up to one (the mixture weights).}
#'  \item{\code{states}}{a vector of length \eqn{N} with entry \eqn{i} being
#'                             the label \eqn{k = 1, \ldots, K} of PLC \eqn{i}}                    
#' } 
#'
#' This is an early release: some function names and arguments might/will
#' (slightly) change in the future, so regularly check with new package updates.
#' @keywords package
#' @references 
#' Shalizi, C. R., R. Haslinger, J.-B. Rouquier, K. L. Klinkner, and C. Moore (2006). 
#' ``Automatic filters for the detection of coherent structure in spatiotemporal systems.''
#' Physical Review E 73, 036104
#' 
#' Shalizi, C. R., K. L. Klinkner, and R. Haslinger (2004a). 
#' ``Quantifying self-organization with optimal predictors.'' 
#' Physical Review Letters 93, 118701.
#' 
#' @seealso 
#' The main functions in this package are
#' \itemize{ 
#'  \item \code{\link{states2LSC}} to estimate LSC from the state space, and
#'  \item \code{\link{LICORS2LSC}} which is a wrapper for estimating LSC from 
#'        a \code{"LICORS"} class estimate.
#' }
#' Since pattern discovery without visualization is only of very limited use, the 
#' \code{\link{plot.LSC}} function shows informative plots for \eqn{(1+1)D} and
#' \eqn{(2+1)D} systems.
#' 
#' @import LICORS RColorBrewer fields gam Matrix
#' 
#' @examples
#' ## known predictive state space with a state-vector
#' data(contCA00)
#' ll <- states2LSC(states = contCA00$predictive_states - min(contCA00$predictive_states) + 1)
#' image2(ll, density = TRUE, legend = FALSE)
#'
#' # An example using estimates from LICORS
#' \dontrun{
#' example(LICORS) # this will give an object 'mod' of class 'LICORS'
#' image2(LICORS2LSC(mod))
#' }
#'
 
NULL
