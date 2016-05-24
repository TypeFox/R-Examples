#' @title Genetic algorithm combined with PLS regression (GA-PLS)
#'
#' @description A subset search algorithm inspired by biological
#' evolution theory and natural selection.
#'
#' @param y vector of response values (\code{numeric} or \code{factor}).
#' @param X numeric predictor \code{matrix}.
#' @param GA.threshold the change for a zero for mutations and initialization (default = 0.5).
#' @param iters the number of iterations (default = 5).
#' @param popSize the population size (default = 100).
#'
#' @details 1. Building an initial population of variable sets by setting bits for each variable
#'  randomly, where bit '1' represents selection of corresponding variable while '0' presents
#'  non-selection. The approximate size of the variable sets must be set in advance.
#'  
#'  2. Fitting a PLSR-model to each variable set and computing the performance by, for instance,
#'  a leave one out cross-validation procedure.
#'  
#'  3. A collection of variable sets with higher performance are selected to survive until the
#'   next "generation".
#'   
#'  4. Crossover and mutation: new variable sets are formed 1) by crossover of selected
#'  variables between the surviving variable sets, and 2) by changing (mutating) the bit
#'  value for each variable by small probability.
#'  
#'  5. The surviving and modified variable sets form the population serving as input to point 2.
#'  
#' @return Returns a vector of variable numbers corresponding to the model 
#' having lowest prediction error.
#'
#' @author Tahir Mehmood, Kristian Hovde Liland, Solve Sæbø.
#'
#' @references K. Hasegawa, Y. Miyashita, K. Funatsu, GA strategy for variable selection
#'  in QSAR studies: GA-based PLS analysis of calcium channel antagonists, Journal of Chemical
#'  Information and Computer Sciences 37 (1997) 306-310.
#'
#' @seealso \code{\link{VIP}} (SR/sMC/LW/RC), \code{\link{filterPLSR}}, \code{\link{spa_pls}}, 
#' \code{\link{stpls}}, \code{\link{truncation}}, \code{\link{bve_pls}}, \code{\link{mcuve_pls}},
#' \code{\link{ipw_pls}}, \code{\link{ga_pls}}, \code{\link{rep_pls}}.
#'
#' @examples
#' data(gasoline, package = "pls")
#' # with( gasoline, ga_pls(octane, NIR) ) # Time consuming
#'
#' @importFrom genalg rbga.bin
#' @export
ga_pls<- function(y,X, GA.threshold=0.5, iters=5, popSize=100){ 
  evalFunc <- evaluateX
  monitorFunc <- monitor
  n <- ncol(X)
  putPLSVarSelEnv("X", X)
  putPLSVarSelEnv("y", y)
  gapls.results <- rbga.bin(size=n, zeroToOneRatio=GA.threshold, evalFunc=evalFunc, monitorFunc=monitorFunc, popSize=popSize, iters=iters)
  bestChro      <- which(gapls.results$population[which.max(rowSums( gapls.results$population))[1],] == 1)
  return (list(ga.selection=bestChro))
}

evaluateX <- function(chromosome=c()) {
  returnVal = 100
  minComp <- 2
  if (sum(chromosome) < minComp) { 
    returnVal
  } else {
    X <- getPLSVarSelEnv("X")
    y <- getPLSVarSelEnv("y")
    Xtr <- X[ ,chromosome == 1];
    pls.model <- plsr(y~Xtr, validation="CV", ncomp=min(5,sum(chromosome)))
    returnVal <- pls.model$val$PRESS[pls.model$val$ncomp-(minComp-1)]
    returnVal
  }
}

monitor <- function(obj) {
  minEval <- min(obj$evaluations);
  filter  <- obj$evaluations == minEval;
  bestObjectCount <- sum(rep(1, obj$popSize)[filter]);
  # Dealing with the situation that more than one object is best
  if (bestObjectCount > 1) {
    bestSolution <- obj$population[filter,][1,];
  } else {
    bestSolution <- obj$population[filter,];
  }
  outputBest <- paste(obj$iter, " #selected=", sum(bestSolution)," Best (Error=", minEval, "): ", sep="")
  for (var in 1:length(bestSolution)) {
    outputBest <- paste(outputBest, bestSolution[var], " ", sep="")
  }
  outputBest <- paste(outputBest, "\n", sep="")
  #	plot(obj, type="hist");
  #	cat(outputBest);
}

