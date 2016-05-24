require(coda)
require(stringr)
require(roxygen2)
#' A package for diagnostics on MCMC chains using Hellinger Distance.
#' 
#' \tabular{ll}{
#' Package: \tab bmk \cr
#' Type: \tab Package \cr
#' Version: \tab 0.1 \cr
#' Date: \tab 2012-09-28 \cr
#' License: \tab GPL (>=2) \cr
#' LazyLoad: \tab yes \cr
#' }
#' 
#' Using Hellinger Distance to determine the distances between MCMC samples for convergence diagnostics
#' and sensitivity studies.
#' 
#' \code{\link{bmksensitive}} compares two identically dimensioned MCMC chains using the Hellinger distance.  
#' The Hellinger distance is calculated for all parameters so one can identify which parameters may be sensitive
#' to changes in prior distribution specification.
#' 
#' \code{\link{bmkconverge}} calculates the Hellinger distance between consecutive bins of MCMC samples.  This can
#' be used to determine when MCMC chains have converged.  Once the Hellinger distance between the consecutive bins of 
#' MCMC samples is less than a specified threshold one can be assured that the samples are similar and hence converged.
#' 
#' \code{\link{bmksummary}} calculates summaries of MCMC samples such as mean, standard deviation, median, 2.5, 97.5 percent 
#' quanties, Gelman-Rubin convergence diagnostic, effective sample size, minimum and maximum Hellinger distances across all chains
#' for each parameter.
#' 
#' @name bmk-package
#' @aliases bmk
#' @docType package
#' @title  MCMC convergence diagnostics using the Hellinger distance
#' @author Edward L. Boone \email{elboone@@vcu.edu} and Matthew Krachey \email{matthewkrachey@yahoo.com}
#' @import coda, plyr, functional
#' @keywords MCMC
#' @references Boone EL, Merrick JR and Krachey MJ.  
#'   A Hellinger distance approach to MCMC diagnostics.
#'   Journal of Statistical Computation and Simulation, 
#'   \code{DOI:10.1080/00949655.2012.729588}.

