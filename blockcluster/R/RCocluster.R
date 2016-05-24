#' @include coclusterStrategy.R
NULL

#' Co-Clustering Package	
#' 
#' This package performs Co-clustering of binary, contingency, continuous and
#' categorical data-sets. 
#' 
#' This package performs Co-clustering of binary, contingency, continuous and
#' categorical data-sets with utility functions to visualize the Co-clustered
#' data. The package contains a set of functions \code{\link{coclusterBinary}},
#' \code{\link{coclusterCategorical}}, \code{\link{coclusterContingency}},
#' \code{\link{coclusterContinuous}} which perform Co-clustering on various
#' kinds of data-sets and return object of appropriate class (refer to
#' documentation of these functions. The package also contains function
#' \code{\link{coclusterStrategy}} (see documentation of function to know
#' various slots) which returns an object of class \code{\linkS4class{strategy}}.
#' This object can be given as input to co-clustering functions to control
#' various Co-clustering parameters. Please refer to testmodels.R file which is
#' included in "test" directory to see examples with various models and
#' simulated data-sets.
#' 
#' The package also provide utility functions like summary() and plot() to
#' summarize results and plot the original and Co-clustered data respectively.
#' 
#' 
#' @examples 
#' 
#' ## Simple example with simulated binary data
#' ## load data
#' data(binarydata)
#' ## usage of coclusterBinary function in its most simplest form
#' out<-coclusterBinary(binarydata,nbcocluster=c(2,3))
#' #" Summarize the output results
#' summary(out)
#' ## Plot the original and Co-clustered data 
#' plot(out)
#' 
#' @import methods
#' 
#' @name blockcluster
#' @rdname blockcluster
#' 
NULL

#' 
#' An EM strategy to obtain a good optimum.
#' 
#' In Co-clustering, there could be many local optimal where the algorithm may
#' get struck resulting in sub-optimum results. Hence we applied a strategy
#' called XEM strategy to run the EM algorithm. The various steps are defined
#' as follows:
#' 
#' \describe{
#' \item{Step-1, "xem" step:}{Do several runs of: "initialization followed by
#' short run of algorithm (few iterations/high tolerance)". This parameter is
#' named as "nbxem" in \code{\link{coclusterStrategy}} function. Default value
#' is 5. We call this step as xem step.}
#' \item{Step-2, "XEM" step:}{Select the best result of step 1 and make long run
#' of Algorithm(high iterations/low tolerance).We call this step as XEM step.}
#' \item{Step-3}{Repeat step 1 and 2 several times and select the best result.
#' The number of repetitions can be modified via parameter "nbtry" of
#' \code{\link{coclusterStrategy}} function. Default value is 2. } 
#' }
#' 
#' @name XEMStrategy
#' @rdname XEMStrategy
#' 
NULL

#' Summary function.
#' 
#' This function gives the summary of output from \code{coclusterBinary},
#' \code{coclusterCategorical}, \code{coclusterContingency},
#' \code{coclusterContinuous}.
#' 
#' @param object output object from \code{\link{coclusterBinary}},
#' \code{\link{coclusterCategorical}}, \code{\link{coclusterContingency}},
#' \code{\link{coclusterContinuous}}.
#' @param ... Additional argument(s) . Currently there is no additional arguments.  
#' @name summary
#' @rdname summary-methods
#' @docType methods
#' @exportMethod summary
#' 
NULL

#' Plot function.
#' 
#' This function plot the original and Co-clustered data-sets.
#' 
#' @param x output object from \code{\link{coclusterBinary}},
#' \code{\link{coclusterCategorical}}, \code{\link{coclusterContingency}},
#' \code{\link{coclusterContinuous}}.
#' @param y Ignored
#' @param ... Additional argument(s). Currently we support two additional argument.
#' "asp": If this is set to TRUE the original aspect ratio is conserved. By
#' default "asp" is FALSE. "type" : This is the type of plot which is either
#' "cocluster" or "distribution". The corresponding plots are Co-clustered data
#' and distributions and mixture densities for Co-clusters respectively.
#' Default is "cocluster" plot.
#' 
#' @importFrom graphics plot
#' @name plot
#' @rdname plot-methods
#' @docType methods
#' @exportMethod plot
NULL

