

#' Allows for network reconstruction and changepoint detection.
#' 
#' This package runs an MCMC simulation to reconstruct networks from time
#' series data, using a non-homogeneous, time-varying dynamic Bayesian network.
#' Networks segments and changepoints are inferred concurrently, and
#' information sharing priors provide a reduction of the inference uncertainty.
#' 
#' \tabular{ll}{ Package: \tab EDISON\cr Type: \tab Package\cr Version: \tab
#' 1.1.1\cr Date: \tab 2016-03-30\cr License: \tab GPL-2\cr LazyLoad: \tab yes\cr
#' }
#' 
#' @name EDISON-package
#' @aliases EDISON-package EDISON
#' @docType package
#' @import MASS
#' @import corpcor
#' @author Frank Dondelinger, Sophie Lebre
#' 
#' Maintainer: Frank Dondelinger <fdondelinger.work@@gmail.com>
#' @seealso \code{\link[corpcor:corpcor-package]{corpcor}}
#' 
#' @references Dondelinger et al. (2012), "Non-homogeneous dynamic Bayesian
#' networks with Bayesian regularization for inferring gene regulatory networks
#' with gradually time-varying structure", Machine Learning.
#' 
#' Husmeier et al. (2010), "Inter-time segment information sharing for
#' non-homogeneous dynamic Bayesian networks", NIPS.
#' @keywords package
#' @examples
#' 
#' # Generate random gene network and simulate data from it
#' dataset = simulateNetwork(l=25)
#' 
#' # Run MCMC simulation to infer networks and changepoint locations
#' result = EDISON.run(dataset$sim_data, num.iter=500)
#' 
#' # Calculate posterior probabilities of changepoints
#' cps = calculateCPProbabilities(result)
#' 
#' # Calculate marginal posterior probabilities of edges in the network
#' network = calculateEdgeProbabilities(result)
#' 
#' 
NULL



