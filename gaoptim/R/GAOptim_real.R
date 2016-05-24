#' Genetic Algorithm setup
#' 
#' Setup a \code{GAReal} object that can be used to perform a real-based optimization.
#' 
#' This is the function used to configure and fine-tune a real-based optimization.
#' The basic usage requires only the \code{FUN} parameter (function to be maximized),
#' together with the \code{lb} and \code{ub} parameters (lower and upper search domain),
#' all the other parameters have sensible defaults.
#' 
#' The parameters \code{selection}, \code{crossover} and \code{mutation} can also take a custom
#' function as argument, which needs to be in the appropriate format (see the
#' examples). The text below explains the default behaviour for these parameters, which will
#' be usefull if you want to override one or more genetic operators.
#' 
#' \itemize{
#' \item \code{selection}: The \code{fitness} option performs a \emph{fitness-proportionate}
#' selection, so that the fittest individuals will have greater chances of being selected.
#' If you choose this option, the value returned by \code{FUN} (the fitness value)
#' should be \strong{non-negative}.
#' The \code{uniform} option will randomly sample the individuals to mate, regardless of
#' their fitness value. See the examples if you want to implement a custom selection function.
#' }
#' 
#' \itemize{
#' \item \code{crossover}: The \code{blend} option will perform a linear combination 
#' of the individuals DNA, effectively introducing new information into the resulting offspring.
#' For details, see \emph{Practical genetic algorithms} in the references. The \code{two.points}
#' option will perform the classic 2-point crossover. See the examples if you need to implement
#' a custom crossover function.
#' }
#' 
#' \itemize{
#' \item \code{mutation}: The default implementation will uniformly sample \code{n} mutation
#' points along the population matrix, where \code{n} is given by \code{mutRate * popSize * nvars} and
#' \code{nvars} is the number of variables in your problem. Each sampled \emph{locus} will be
#' replaced by a random-uniform number between 0 and 1. See the examples to learn how to use
#' a custom mutation function.
#' }
#' 
#' @param FUN The fitness function, which should take a vector as argument and return a numeric
#' value (See details).
#' @param lb A numeric vector specifying the lower bounds for the search domain.
#' @param ub A numeric vector specifying the upper bounds for the search domain.
#' @param popSize The population size.
#' @param mutRate The mutation rate, a numeric value between 0 and 1. When implementing a
#' custom mutation function, this value should be one of the parameters (see details and examples).
#' @param cxRate The crossover rate, a numeric value between 0 and 1. This parameter specifies
#' the probability of two individuals effectively exchange DNA during crossover. In case the
#' individuals didn't crossover, the offspring is a exact copy of the parents. When implementing
#' a custom crossover function, this value should be one of the arguments (see details and examples).
#' @param eliteRate A numeric value between 0 and 1. The \code{eliteRate * popSize} best-fitted
#' individuals will automatically be selected for the next generation.
#' @param selection The selection operator to be used. You can also implement a custom selection
#' function (see details and examples).
#' @param crossover The crossover operator to be used. You can also implement a custom crossover
#' function (see details and examples).
#' @param mutation The mutation operator to be used. You can also implement a custom mutation
#' function (see details and examples).
#' @return An object of class \code{GAReal}, which you can pass as an argument to \code{plot} or
#' \code{summary}. This object is a list with the following accessor functions:
#' \tabular{ll}{
#' \code{bestFit}: \tab Returns a vector with the best fitness achieved in each generation.\cr
#' \code{meanFit}: \tab Returns a vector with the mean fitness achieved in each generation.\cr
#' \code{bestIndividual}: \tab Returns a vector with the best solution found.\cr
#' \code{evolve(h)}: \tab This is the function you call to evolve your population.
#' \cr \tab You also need to specify the number of generations to evolve.\cr
#' \code{population}: \tab Returns the current population matrix.
#' }
#' @export
#' @examples
#' # Maximize a trivial 5 variable function
#' # The function and search-space below will be used for all examples
#' 
#' fitness.FUN = function(x) sum(x)
#' lb = c(0, 0, 0, 0, 0)
#' ub = c(10, 10, 10, 10, 10)
#' 
#' ga1 = GAReal(fitness.FUN, lb, ub)
#' ga1$evolve(200)
#' plot(ga1)
#' 
#' # A custom selection example
#' selec.FUN = function(population, fitnessVec, nleft)
#' {
#'  # population - The population matrix
#'  # fitnessVec - The corresponding fitness vector for the population matrix
#'  # nleft - The number of individuals you should select
#'  
#'  half = as.integer(nleft/2)
#'  remain = nleft - half
#'  idxs = 1:nrow(population)
#'  
#'  # pick half using fitness-proportionate
#'  rowIdxs = sample(idxs, half, replace = TRUE, prob = fitnessVec)
#'  # pick the other half randomly
#'  rowIdxs = c(rowIdxs, sample(idxs, remain, replace = TRUE))
#'  
#'  # Just return the nLeft selected row indexes
#'  return(rowIdxs)
#' }
#' 
#' ga2 = GAReal(fitness.FUN, lb, ub, selection = selec.FUN)
#' ga2$evolve(200)
#' summary(ga2)
#' 
#' # A custom crossover example
#' crossover.FUN = function(parent1, parent2, prob)
#' {
#'  # parent1, parent2 - The individuals to crossover
#'  # prob - The probability of a crossover happen (cxRate parameter)
#'  
#'  # Respect the cxRate parameter: if DNA is not exchanged, just return the parents
#'  if (runif(1) > prob)
#'    return(matrix(c(parent1, parent2), nrow = 2, byrow = TRUE))
#'    
#'  # A simple uniform crossover - just swap the 'genes' with a probability of 0.5
#'  for (i in 1:length(parent1))
#'  {
#'    if (runif(1) > 0.5)
#'    {
#'      tempval = parent1[i]
#'      parent1[i] = parent2[i]
#'      parent2[i] = tempval      
#'    }
#'  }
#'  # You should return a matrix in this format
#'  return(matrix(c(parent1, parent2), nrow = 2, byrow = TRUE))
#' }
#' 
#' ga3 = GAReal(fitness.FUN, lb, ub, crossover = crossover.FUN)
#' ga3$evolve(200)
#' plot(ga3)
#'
#' # A custom mutation example
#' mutation.FUN = function(population, nMut)
#' {
#'  # population - The population matrix to apply mutation
#'  # nMut - The number of mutations you supposed to apply, according to mutRate
#'  
#'  rows = sample(1:nrow(population), nMut, replace = TRUE)
#'  cols = sample(1:ncol(population), nMut, replace = TRUE)
#'  noise = (runif(nMut))^2
#'  
#'  # extract the matrix indexes
#'  ext = matrix(c(rows, cols), nMut, 2)
#'  population[ext] = noise
#'  return(population)
#' }
#' 
#' ga4 = GAReal(fitness.FUN, lb, ub, mutation = mutation.FUN) 
#' ga4$evolve(200)
#' summary(ga4)
#'  
#' @references Randy L. Haupt, Sue Ellen Haupt (2004). Practical genetic 
#' algorithms - 2nd ed.
#' @references Michalewicz, Zbigniew. Genetic Algorithms + Data Structures = Evolution
#' Programs - 3rd ed.
#' 
GAReal = function (FUN, lb, ub, popSize = 100, mutRate = 0.01, cxRate = 0.9, eliteRate = 0.4,
                     selection = c('fitness', 'uniform'), crossover = c('blend', 'two.points'),
                     mutation = c('noise'))
{		
  
  # Basic arg exception check #
  if (! is.numeric(popSize) || popSize < 4)
    stop("Please set 'popSize' to an integer value greater than 3.")
  
  if (! is.numeric(mutRate) || mutRate < 0 || mutRate > 1)
    stop("Please set 'mutRate' to a value in the range [0, 1].")
  
  if (! is.numeric(cxRate) || cxRate < 0 || cxRate > 1)
    stop("Please set 'cxRate' to a value in the range [0, 1].")
  
  if (! is.numeric(eliteRate) || eliteRate < 0 || eliteRate >= 1)
    stop("Please set 'eliteRate' to a value in the range [0, 1[.")
  # Basic arg exception check #
  
  currentPopulation = NULL		
  bestFitnessVec = numeric()
  meanFitnessVec = numeric()
  elite = max(0, 2 * as.integer(eliteRate * popSize * 0.5))
  popSize = 2 * as.integer(popSize * 0.5)
  nvars = length(lb)
  bestCX = rep(0, nvars)
  bestFit = NULL
  mutations = as.integer(mutRate * popSize * nvars)
  iter = 0
  newPopulation = matrix(0, nrow = popSize, ncol = nvars)
  
  ############### BEG selection function definitions #########################################
  selection.FUN = NULL
  if (is.function(selection))
  	selection.FUN = selection
  else
  	selection.type = switch(match.arg(selection), fitness = 'fitness', uniform = 'uniform')
  ############### END selection function definitions #########################################
  	
  ############### BEG crossover function definitions #########################################
  twoPointsCrossover = function(x1, x2, prob)
  {	
  	if (runif(1) > prob)
    	return(matrix(c(x1, x2), nrow = 2, byrow = TRUE))

	p12 = sample(1:length(x1), 2)
	idxs = seq(p12[1], p12[2])
	temp1 = x1[idxs]
	x1[idxs] = x2[idxs]
	x2[idxs] = temp1
	matrix(c(x1, x2), nrow = 2, byrow = TRUE)
  }
  
  blendCrossover = function(cr1, cr2, prob)
  {								
    if (runif(1) > prob)                  
    	return(matrix(c(cr1, cr2), nrow = 2, byrow = TRUE))
    
    beta = 0.5
    n = length(cr1)
    i = sample(1:n, 1)
    
    pm = cr1[i]
    pd = cr2[i]			
    ch1 = cr1
    ch2 = cr2			
    ch1[i:n] = cr2[i:n]
    ch2[i:n] = cr1[i:n]
    
    ch1[i] = pm - beta*(pm - pd)
    ch2[i] = pd + beta*(pm - pd)
    
    matrix(c(ch1, ch2), nrow = 2, byrow = TRUE)			
  }

  applyCrossover = function(rowIdxs, M, FUN)
  {
  	FUN.vec = function(rowVector, mat)
  	{
  		FUN(mat[rowVector[1], ], mat[rowVector[2], ], cxRate)
  	}
  	
    m1 = apply(rowIdxs, 1, FUN.vec, mat = M)
    matrix(t(m1), byrow = F, ncol = ncol(M))
  }
    
  crossover.FUN = NULL
  if (is.function(crossover))
  	crossover.FUN = crossover
  else
  	crossover.FUN = switch(match.arg(crossover), blend = blendCrossover, two.points = twoPointsCrossover) 
  ############### END crossover function definitions #########################################
  
  ############### BEG mutation function definitions #########################################
  mutateNoise = function(x, mutations)
  {			
    rows = sample(1:nrow(x), mutations, replace = TRUE)
    cols = sample(1:ncol(x), mutations, replace = TRUE)
    noise = runif(mutations)
    
    ext = matrix(c(rows, cols), mutations, 2)
    x[ext] = noise
    x
  }
  
  	mutation.FUN = NULL
  	if (is.function(mutation))
  		mutation.FUN = mutation
  	else
  		mutation.FUN = switch(match.arg(mutation), noise = mutateNoise)
  ############### END mutation function definitions #########################################
      
  decode = function(x, lb, ub)
  {
    n = nrow(x)
    x * rep(ub - lb, each = n) + rep(lb, each=n) 	
  }
  
  initPopulation = function()
  {		
    if (is.null (currentPopulation))
    {
      eps = 10E-6
      if (length(lb) != length(ub))
        stop('Domain vectors must have the same length.')
      if (any (is.na (lb)) || any( is.na (ub)))
        stop('Missing values not allowed in Domain vectors.')
      if ( any(ub - lb < eps) )
        stop('Values are too close in Domain vectors. Please specify a wider search space')
      
      n = popSize * nvars
      currentPopulation <<- matrix (runif (n), nrow = popSize)				
    }			
  }
  
  initPopulation()
  
  do.evolve = function()
  {
    iter <<- iter + 1      
    decodedPop = decode(currentPopulation, lb, ub)
    fitnessVec = apply(decodedPop, 1, FUN)
    this.best = max(fitnessVec)			
    bestFitnessVec[iter] <<- this.best			
    meanFitnessVec[iter] <<- mean(fitnessVec)
    
    if (is.null(bestFit) || (this.best > bestFit))
    {
      bestFit <<- this.best
      # [1]: pode haver mais de 1 instancia do melhor individuo
      bestCX <<- decodedPop[which(fitnessVec == this.best)[1], ]
    }
    
    nLeft = popSize
    
    if (elite > 0)
    {
      # Maximization problem												
      nLeft = popSize - elite
      newPopulation[1:elite, ] = currentPopulation[order(fitnessVec, decreasing = TRUE)[1:elite], ] 
    }
    
    # crossover selection
    if (! is.null(selection.FUN))
    	popIdxs = selection.FUN(decodedPop, fitnessVec, nLeft)[1:nLeft]
    else
    { 
    	if (identical(selection.type, 'fitness'))
      		probVec = fitnessVec
    	else if (identical(selection.type, 'uniform'))
      		probVec = NULL
      		    	      	
      	popIdxs = sample(1:popSize, nLeft, replace = TRUE, prob = probVec)
    }
    
    popIdxsM = matrix(popIdxs, ncol = 2, byrow = TRUE)     
    offspring = applyCrossover(popIdxsM, currentPopulation, crossover.FUN)
    newPopulation[(elite+1):popSize, ] = offspring
    
    currentPopulation <<- mutation.FUN(newPopulation, mutations)
  }		
  
  objs = list (		
    population = function()
    {
      decode(currentPopulation, lb, ub)
    },								
    
    bestFit = function()
    {
      bestFitnessVec
    },
    
    meanFit = function()
    {
      meanFitnessVec
    },
    
    bestIndividual = function()
    {
      bestCX
    },
    
    evolve = function(h)
    {        	
      if (missing(h))
        stop('Please specify the number of generations to evolve.\n')
      
      length(bestFitnessVec) = length(bestFitnessVec) + h
      length(meanFitnessVec) = length(meanFitnessVec) + h
      invisible(replicate(h, do.evolve()))
    }
    
  )
  
  class(objs) = 'GAReal'
  objs
}

#' Genetic Algorithm plot
#' 
#' A quick way to visualize the GA results.
#' 
#' @param x An object of class \code{GAReal}.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @param main The plot title.
#' @param bestcol The color for the best fitness evolution line
#' @param meancol The color for the mean fitness evolution line
#' @param lwd The line width.
#' @param legend.pos The legend position, as a character vector.
#' @param ... Other parameters (will be ignored).
#' @aliases plot.GAReal
#' @export
#' @examples
#' 
#' ga = GAReal(function(x) sum(x), rep(0, 5), rep(10, 5))
#' ga$evolve(200)
#' plot(ga)
#' 

#' @method plot GAReal
#' @S3method plot GAReal
#' @rdname plot_real
plot.GAReal = function(x, xlab = 'Generation', ylab = 'Fitness', main = 'GA optimization',
		                  bestcol = 'steelblue', meancol = 'tomato', lwd = 2, 
		                  legend.pos = c('bottomright', 'bottom', 'bottomleft',
                      'left', 'topleft', 'top', 'topright', 'right', 'center'), ...)
{
	ymean = x$meanFit()
	if (length(ymean) == 0)
	{
		print(summary(x))
		return(NULL)
	}
	
	ybest = x$bestFit()
	ylim = c(min(ymean, ybest), max(ymean, ybest))
	plot(ybest, col = bestcol, panel.first = grid(col = '#A9A9A9'), xlab = xlab,
	ylab = ylab, main = main, lwd = lwd, type = 'l', ylim = ylim)
	
	lines(ymean, col = meancol, lwd = lwd)
	legend(legend.pos, legend = c('best', 'mean'), col = c(bestcol, meancol), lwd = lwd)
}

#' Genetic Algorithm summary
#' 
#' A simple summary of the GA results.
#' 
#' @return An object of class \code{summaryGAReal}, which is a list that can be inspected or
#' printed on-screen.
#' @param object An object of class \code{GAReal}.
#' @param ... Other parameters (will be ignored).
#' @export
#' @aliases summary.GAReal

#' @method summary GAReal
#' @S3method summary GAReal
#' @rdname summary_real
summary.GAReal = function(object, ...)
{
	n = length(object$bestFit())		
	if (n == 0)
	{
		sm.obj = list(evolved = FALSE)
		class(sm.obj) = 'summaryGAReal'
		return(sm.obj)		
	}
	sm.obj = list(n = n, sm.mean = summary(object$meanFit()),
				sm.best = summary(object$bestFit()),
				best.cx = object$bestIndividual(),
				best.fit = max(object$bestFit()), 
				evolved = TRUE)
	class(sm.obj) = 'summaryGAReal'
	sm.obj
}

#' Print GA results
#' 
#' Prints the GA results.
#' @param x An object of class \code{GAReal} or \code{summaryGAReal}
#' @param ... Other parameters (will be ignored).
#' @export
#' @aliases print.GAReal

#' @method print GAReal
#' @S3method print GAReal
#' @rdname print_real
print.GAReal = function(x, ...)
{
	print(summary(x))
}

#' @method print summaryGAReal
#' @S3method print summaryGAReal
#' @rdname print_real
print.summaryGAReal = function(x, ...)
{
	if (! x$evolved)
	{
		cat('Population ready to evolve.')
		cat('\nPlease, call myGA$evolve(h) to generate results.\n')
	}
	else
	{
		cat('Results for', x$n, 'Generations:')
		cat('\nMean Fitness:\n')
		print(x$sm.mean)
		cat('\nBest Fitness:\n')
		print(x$sm.best)
		cat('\nBest individual:\n')
		print(x$best.cx)
		cat('\nBest fitness value:\n')
		print(x$best.fit)
	}
}