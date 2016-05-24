#' Strategy function
#' 
#' This function is used to set all the parameters for Co-clustering. It returns
#' an object of class  \code{\linkS4class{strategy}} which can be given as input
#' to \code{\link{coclusterBinary}},
#' \code{\link{coclusterCategorical}}, \code{\link{coclusterContingency}},
#' \code{\link{coclusterContinuous}} function. 
#' 
#' @param algo  The valid values for this parameter are "BEM" (Default), "BCEM"
#' and "BSEM".
#' @param stopcriteria It specifies the stopping criteria. It can be based on
#' either relative change in parameters (preffered due to computation reasons)
#' value or relative change in pseudo log-likelihood. Valid criterion values are
#' "Parameter" and "Likelihood". Default criteria is "Parameter".
#' @param initmethod Method to initialize model parameters. The valid values are
#' "cemInitStep", "fuzzyCemInitStep" and "randomInit". For now only one kind of
#' initialization exists for every model currently available in the package.
#' Hence default value for initialization is set according to the model.
#' @param nbinititerations Number of Global iterations used in initialization
#' step. Default value is 10.
#' @param initepsilon Tolerance value used while initialization. Default value
#' is 1e-2.
#' @param nbiterations_int Number of iterations for internal E step. Default value is 5.
#' @param epsilon_int Tolerance value for relative change in Parameter/likelihood
#' for internal E-step. Default value is 1e-2.
#' @param nbtry Number of tries (XEM steps). Default value is 2.
#' @param nbxem Number of xem steps. Default value is 5.
#' @param nbiterationsxem Number of EM iterations used during xem step.
#' Default value is 50.
#' @param nbiterationsXEM Number of EM iterations used during XEM step.
#' Default value is 500.
#' @param epsilonxem Tolerance value used during xem step. Default value is 1e-4.
#' @param epsilonXEM Tolerance value used during XEM step. Default value is 1e-10
#' @param hyperparam Hyper-parameters ("a" and "b") in case of Bayesian settings (only valid
#' for binary and categorical models). Default values are c(1,1) (no prior).
#' 
#' @return Object of class \code{\linkS4class{strategy}}
#' @export
#' 
#' @examples
#' #Default strategy values
#' 
#' strategy<-coclusterStrategy()
#' summary(strategy)
#' 
coclusterStrategy<-function(algo = "BEM",initmethod=character(),stopcriteria = "Parameter",  
		nbiterationsxem = 50, nbiterationsXEM = 500,nbinititerations = 10, initepsilon = 1e-2, nbiterations_int = 5,
		epsilon_int = 1e-2, epsilonxem = 1e-4,epsilonXEM =1e-10, nbtry = 2, nbxem = 5, hyperparam = c(1,1))
{ 
	#create and return object of class strategy
	new("strategy",algo = algo,initmethod = initmethod,stopcriteria = stopcriteria, 
			nbinititerations = nbinititerations,initepsilon = initepsilon, nbiterations_int = nbiterations_int, 
			epsilonxem = epsilonxem,epsilonXEM = epsilonXEM, epsilon_int = epsilon_int, 
			nbtry = nbtry,nbxem = nbxem,nbiterationsxem = nbiterationsxem, nbiterationsXEM = nbiterationsXEM,
      hyperparam=hyperparam)
}

#' strategy class
#' 
#' This class contains all the input parameters to run coclustering.
#' 
#' \describe{
#' \item{algo: }{Algorithm to be use for co-clustering.}
#' \item{stopcriteria: }{Stopping criteria used to stop the algorithm.}
#' \item{initmethod: }{Method to initialize model parameters.}
#' \item{nbinititerations: }{Number of global iterations while running initialization.}
#' \item{initepsilon: }{Tolerance value used while initialization.}
#' \item{nbiterations_int: }{Number of iterations for internal E-step.}
#' \item{epsilon_int: }{Tolerance value for internal E-step.}
#' \item{nbtry: }{Number of tries.}
#' \item{nbxem: }{Number of xem iterations.}
#' \item{nbiterationsxem: }{Number of EM iterations used during xem.}
#' \item{nbiterationsXEM: }{Number of EM iterations used during XEM.}
#' \item{epsilonxem: }{Tolerance value used during xem.}
#' \item{epsilonXEM: }{Tolerance value used during XEM.}
#' \item{hyperparam: }{Hyper-parameters ("a" and "b") in case of Bayesian settings.}
#' }
#' 
#' @rdname coclusterStrategy
#' @name strategy-class
#' @exportClass strategy
#' 

setClass(
		Class="strategy",
		representation = representation(
				algo = "character",
				initmethod = "character",
				stopcriteria = "character",
				nbinititerations = "numeric",
				initepsilon = "numeric",
				nbiterations_int = "numeric",
				epsilon_int = "numeric",
				nbtry = "numeric",
				nbxem = "numeric",
				nbiterationsxem = "numeric",
				nbiterationsXEM = "numeric",
				epsilonxem = "numeric",
				epsilonXEM = "numeric",
				hyperparam = "numeric"
		),
		prototype = prototype(
				algo = character(0),
				stopcriteria = character(0),
				initmethod = character(0),
				nbinititerations = integer(0),
				initepsilon = numeric(0),
				nbiterations_int = integer(0),
				epsilon_int = numeric(0),
				nbtry = integer(0),
				nbxem = integer(0),
				nbiterationsxem = integer(0),
				nbiterationsXEM = integer(0),
				epsilonxem = numeric(0),
				epsilonXEM = numeric(0),
				hyperparam = numeric(0)
		)
)
			

#' @rdname getter-methods
#' @aliases [,strategy-method
setMethod(
		f = "[",
		signature = "strategy",
		definition = function(x,i,j,drop) {
			switch(EXPR=i,
					"stopcriteria"={return (x@stopcriteria)}, 
					"initmethod"={return (x@initmethod)},
					"nbcocluster"={return (x@nbcocluster)}, 
					"nbinititerations"={return (x@nbinititerations)},
					"initepsilon" = {return (x@initepsilon)},
					"nbiterations_int"={return (x@nbiterations_int)},  
					"epsilon_int"={return (x@epsilon_int)},  
					"nbtry"={return (x@nbtry)},
					"nbxem"={return (x@nbxem)},
					"nbiterationsxem"={return (x@nbiterationsxem)},
					"nbiterationsXEM"={return (x@nbiterationsXEM)},
					"epsilonxem"={return (x@epsilonxem)},
					"epsilonXEM"={return (x@epsilonXEM)},
          "hyperparam"={return (x@hyperparam)},
					stop("Invalid slot name.")
			)			
		}
)

#' @rdname summary-methods
#' @aliases summary summary,strategy-method
setMethod(
    f="summary",
    signature = "strategy",
    definition = function(object,...) {
      cat("******************************************************************\n")
      cat("Algorithm: ",object@algo)
      cat("\nInitialization method(There is no default value): ",object@initmethod)
      cat("\nStopping Criteria: ",object@stopcriteria)
      cat("\n\nVarious Iterations")
      cat("\n******************")
      cat("\nNumber of global iterations while running initialization: ",object@nbinititerations)
      cat("\nNumber of iterations for internal E-step: ",object@nbiterations_int)
      cat("\nNumber of EM iterations used during xem: ",object@nbiterationsxem)
      cat("\nNumber of EM iterations used during XEM: ",object@nbiterationsXEM)
      cat("\nNumber of xem iterations: ",object@nbxem)
      cat("\nNumber of tries: ",object@nbtry)
      cat("\n\nVarious epsilons")
      cat("\n****************")
      cat("\nTolerance value used while initialization: ",object@initepsilon)
      cat("\nTolerance value for internal E-step: ",object@epsilon_int)
      cat("\nTolerance value used during xem: ",object@epsilonxem)
      cat("\nTolerance value used during XEM: ",object@epsilonXEM)
      cat("\nHyper-parameters: ",object@hyperparam)
      cat("\n******************************************************************\n")
    }
)
