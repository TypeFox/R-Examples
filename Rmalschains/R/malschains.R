#' This is a function that initializes and sets the parameters of the algorithm. 
#' It generates a list of parameters, to be used with the \code{\link{malschains}} function. 
#'
#' @title Sets and initializes the main parameters of the algorithm
#' @param popsize The population size of the evolutionary algorithm.
#' @param ls The local search method. Currently implementend are \code{cmaes}, \code{sw}, \code{simplex}, and \code{ssw}.
#' Usually, the \code{cmaes} local search strategy will give good results. However, it does not scale well with the problem size. 
#' So, if performance is needed, the \code{sw} strategy is a better choice. If the problem is high-dimensional, the \code{ssw} strategy is promising,
#' which selects randomly 20\% of the variables for optimization.
#' @param istep The number of iterations of the local search. I.e., if the local search is started or re-started on an individual, it will be executed for an 
#' \code{istep} number of iterations. 
#'  This parameter depends on the local search used. For \code{cmaes}, usually an \code{istep} of 300 is
#' a good choice. For the other local search methods, an \code{istep} of 100 performs better.
#' @param effort A value between 0 and 1 which gives the ratio between the amount of evaluations that are used for the local search and
#' for the evolutionary algorithm, respectively. A higher effort means more evaluations for the evolutionary algorithm. So, if exploration 
#' of the search space is more important than finding local optima, \code{effort} is to be chosen higher.
#' @param alpha The alpha parameter from crossover BLX-alpha. A lower value (< 0.3) reduces diversity, a higher value increases diversity.
#' @param optimum The optimum to achieve. The default is zero, as in many minimization problems a value of zero can be considered optimal. 
# @param threshold A threshold which defines how far away from the optimum the solution will still be treated as optimal.
#' @param threshold A threshold which defines for the local search how much improvement is considered as no improvement. If this value is chosen
#' too low (zero), then the local search will usually always try to improve on the best individual, even if it is already located very close to a local optimum. 
#' @param lsOnly Apply only the local search algorithm, and not MA-LS-Chains
#' @param lsParam1 First (optional) parameter. Currently, if local search is cmaes, this is the parameter popsize/lambda of cmaes. If it is not set or not positive, cmaes will calculate this automatically using a heuristic.
#' @param lsParam2 Second (optional) parameter. Currently, if local search is cmaes, this is the parameter parentssize/mu of cmaes. If it is not set, not positive, or not smaller lambda, cmaes will calculate this automatically using a heuristic.
#' @references 
#' 
#' Molina, D., Lozano, M., Sánchez, A.M., Herrera, F.
#' Memetic algorithms based on local search chains for large scale continuous optimisation problems: MA-SSW-Chains
#' (2011) Soft Computing, 15 (11), pp. 2201-2220. 
#'
#' Molina, D., Lozano, M., García-Martínez, C., Herrera, F.
#' Memetic algorithms for continuous optimisation based on local search chains
#' (2010) Evolutionary Computation, 18 (1), pp. 27-63.
#' 
#' @export
malschains.control <- function(popsize=50, ls="cmaes", istep=500, effort=0.5, alpha=0.5, optimum=-Inf, threshold=1e-8, lsOnly=FALSE, lsParam1=0, lsParam2=0) {
  
  if(popsize %% 10 != 0) {
    warning("Only population sizes divisible by 10 are supported. Using rounded value: ", popsize)
    popsize <- round(popsize / 10, digits=0) * 10
  }
  if(popsize == 0) popsize <- 10
  
  list(popsize=popsize, ls=ls, istep=istep, effort=effort, alpha=alpha, optimum=optimum, threshold=threshold, lsOnly=lsOnly, lsParam1=lsParam1, lsParam2=lsParam2)
}


#' This is the main function of the package. It minimizes the output of the function fn (for maximization, change the sign of the output of fn). 
#' 
#' The output of the function when run with \code{verbosity=2} is the following:
#' 
#' \itemize{
#' \item \code{EA::PopFitness} The fitness of the best, the one at the 1st quartile, the one at the 3rd quartile, and the worst individual.
#' \item \code{EA::Improvement} Improvement of the individuals at the according ranked positions in the population (best, 1st quartile, 3rd quartile, worst).
#' \item \code{LS} The number of the individual which is improved on (in braces), its fitness before and after application of the LS procedure, and their difference. 
#' \item \code{EABest} If the best fitness present in the population changed: same as \code{LS}.
#' }
#' 
#' @title Perform optimization with the MA-LS-Chains algorithm
#' @param fn The function to minimize.
#' @param lower The lower bound (or bounds) of the search domain. 
#' @param upper The upper bound (or bounds) of the search domain.
#' @param dim The dimension of the problem (if \code{lower} and \code{upper} are vectors it is not needed).
#' @param maxEvals The maximal number of evaluations of the fitness function. 
#' @param verbosity Set the verbosity level. Currently, meaningful values are 0, 1, 2
#' @param initialpop An initial population for the evolutionary algorithm can be submitted (as a matrix). Here, prior knowledge
#' can be introduced to get better results from the algorithm.
#' @param control A list containing the main options of the algorithm. See \code{\link{malschains.control}}.
#' @param env The environment in which to evaluate the fitness function. If not given, it is generated.
#' @param seed A seed value for the random number generator.
#' @return the function returns a list containing the best individual, \code{sol}, and its \code{fitness}. 
#' Furthermore, it contains some information on the optimization process, which can be seen using \code{\link{print.malschains}}.
#' @references 
#' 
#' Molina, D., Lozano, M., Sánchez, A.M., Herrera, F.
#' Memetic algorithms based on local search chains for large scale continuous optimisation problems: MA-SSW-Chains
#' (2011) Soft Computing, 15 (11), pp. 2201-2220. 
#'
#' Molina, D., Lozano, M., García-Martínez, C., Herrera, F.
#' Memetic algorithms for continuous optimisation based on local search chains
#' (2010) Evolutionary Computation, 18 (1), pp. 27-63.
#' 
#' @export
malschains <- function(fn, lower, upper, dim, maxEvals=10*control$istep, verbosity=2, initialpop = NULL, control=malschains.control(), seed=NULL, env) {
  
  dimv = length(lower)
  stopifnot(length(lower)==length(upper))
  
  if (missing(dim)) {
    dim = dimv 
    stopifnot(dim >= 1)
  }
  else {
    stopifnot(dim >= 1)
    lower = rep(lower, dim)
    upper = rep(upper, dim)
  }
  
  istep <- round(control$istep)
  effort <- as.numeric(control$effort)
  alpha <- as.numeric(control$alpha)
  threshold <- as.numeric(control$threshold)
  optimum <- as.numeric(control$optimum)
  lsOnly <- control$lsOnly
  lsParam1 <- as.numeric(control$lsParam1)
  lsParam2 <- as.numeric(control$lsParam2)
  
  # @param minimize boolean (TRUE indicates that it should minimize, FALSE in other case)
  minimize <- TRUE
  
  if(!is.null(initialpop)) {
    if(!is.matrix(initialpop))
      initialpop <- matrix(initialpop, ncol=dim)
  }
  
  if(missing(env)) env <- new.env()
  
  if((maxEvals < istep) & !lsOnly) {
    warning("maxEvals cannot be smaller that istep, setting to istep") 
    maxEvals <- istep
  }
  
  #use the cmaes version which respects the seed
  if (control$ls == "cmaes") control$ls <- "cmaesmyrandom"
  if(is.null(seed)) seed <- 1e8*runif(1)
  
  res <- .Call( "RmalschainsWrapper", fn, dim, lower, upper, env, control$popsize, maxEvals, control$ls, verbosity, istep, effort, 
        alpha, optimum, threshold, minimize, initialpop, seed, lsOnly, lsParam1, lsParam2, PACKAGE = "Rmalschains" )
    
  class(res) <- "malschains" 
  res  
}

#' Print out some characteristics of a \code{\link{malschains}} result.
#' The result shows besides the best solution and its fitness the total number of evaluations spent for both EA and LS, 
#' the ratio of the spent evaluations (also called effort), the ratio of total improvement of the fitness,
#' the percentage of times that application of the EA/LS yielded improvement, and some timing results in milliseconds. 
#'  
#' @title Generic print function for malschains results
#' @param x the \code{\link{malschains}} result
#' @param ... additional function parameters (currently not used)
#' @export
#' @S3method print malschains
#' @method print malschains
# @rdname malschains
print.malschains <- function(x, ...) {
  if(!inherits(x, "malschains")) stop("not a legitimate malschains result")
  
  cat(sprintf("NumTotalEvalEA: %d\n", x$numEvalEA))
  cat(sprintf("NumTotalEvalLS: %d\n", x$numEvalLS))  
  
  ratio_effort <- x$numEvalEA/(x$numEvalEA+x$numEvalLS)
  cat(sprintf("RatioEffort EA/LS: [%.0f/%.0f]\n", 100*ratio_effort, 100*(1-ratio_effort)))
  
  ratio_alg <- x$improvementEA/(x$improvementEA+x$improvementLS)
  cat(sprintf("RatioImprovement EA/LS: [%.0f/%.0f]\n", 100*ratio_alg, 100*(1-ratio_alg)))
  #cat(sprintf(("Restarts: %d\n", restarts))
  
  if((x$numTotalEA != 0) && (x$numTotalLS != 0)) {
    cat(sprintf("PercentageNumImprovement[EA]: %d%%\n", round((x$numImprovementEA*100)/x$numTotalEA)))
    cat(sprintf("PercentageNumImprovement[LS]: %d%%\n", round((x$numImprovementLS*100)/x$numTotalLS)))    
  }
  
  cat(sprintf("Time[EA]: %.2f\n", x$timeMsEA))
  cat(sprintf("Time[LS]: %.2f\n", x$timeMsLS))
  cat(sprintf("Time[MA]: %.2f\n", x$timeMsMA))
  cat(sprintf("RatioTime[EA/MA]: %.2f\n", 100*x$timeMsEA/x$timeMsMA))
  cat(sprintf("RatioTime[LS/MA]: %.2f\n", 100*x$timeMsLS/x$timeMsMA))
  
  
  cat("Fitness:\n",sep="")
  print(x$fitness)
  cat("Solution:\n",sep="") 
  print(x$sol)
  
  invisible(x)
}
