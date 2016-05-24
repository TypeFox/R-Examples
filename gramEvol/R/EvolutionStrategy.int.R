EvolutionStrategy.int <- function(genomeLen, codonMin, codonMax,
                                  genomeMin=rep.int(codonMin, genomeLen), 
                                  genomeMax=rep.int(codonMax, genomeLen),
                                  suggestion=NULL,
                                  popSize=4, newPerGen = 4,
                                  iterations=500, terminationCost=NA,
                                  mutationChance=1/(genomeLen+1),
                                  monitorFunc=NULL, 
                                  evalFunc,
                                  allowrepeat = TRUE,
                                  showSettings=FALSE, verbose=FALSE,
                                  plapply = lapply) {
  # Optimizes an Integer chromosome using Evolutionary strategy
  #
  # popSize          = the offspring population size
  # iterations       = number of generations
  # terminationCost  = The cost (error) that if reached, the GA should termiante
  # mutationChance   = chance that a var in the string gets mutated
  #
  is.verbose = verbose
  verbose = function(...) { if (is.verbose) cat(...)}
  
  if (is.null(evalFunc)) {
    stop("A evaluation function must be provided. See the evalFunc parameter.");
  }
  
  stopifnot(genomeLen > 1)
  
  # do a variaty of sanity checks first
  verbose("Testing the sanity of parameters...\n");
  if (length(genomeMin) != length(genomeMax)) {
    stop("The vectors genomeMin and genomeMax must be of equal length.");
  }
  if (iterations < 1) {
    stop("The number of iterations must be at least 1.");
  }
  if ((mutationChance < 0) | (mutationChance  > 1)) {
    stop("mutationChance must be between 0 and 1.");
  }
  
  if (showSettings) {
    verbose("The start conditions:\n");
    result = list(genomeMin=genomeMin, genomeMax=genomeMax, suggestions=suggestion,
                  popSize=popSize, iterations=iterations,
                  mutationChance=mutationChance);
    class(result) = "rbga";
    
    cat(summary(result));
  } else {
    verbose("Not showing GA settings...\n");
  }
  
  ##########
  # Creation
  if (!is.null(suggestion)) {
    verbose("Adding suggestions to first population...\n");
    suggestionCount = 1
    parent = suggestion
  } else {
    verbose("Starting with random values in the given domains...\n");
    suggestionCount = 0
    
    parent = ga.new.chromosome(genomeLen, genomeMin, genomeMax, allowrepeat)
    parentEval = NA
  }
  
  ############################################################################
  # do iterations
  bestEvals = rep(NA, iterations);
  meanEvals = rep(NA, iterations);
  
  totalPopulation = 1 + popSize + newPerGen
  
  for (iter in 1:iterations) {
    verbose(paste("Starting iteration", iter, "\n"));
    
    ############################################################################
    # Copy the parent to popultion
    population = matrix(rep(parent, totalPopulation), nrow=totalPopulation, byrow = TRUE)
    evalVals = rep(NA, totalPopulation);
    evalVals[1] = parentEval
    
    ############################################################################
    # Mutation
    if (mutationChance > 0 & popSize > 0) {
      verbose("  applying mutations... ");
      mutationCount = 0;
      for (object in 2:popSize) { # don't mutate the parent
        
        dempeningFactor = 1#(iterations-iter)/iterations

        mutResult <- ga.mutation(population[object,], mutationChance, genomeLen, 
                                 genomeMin, genomeMax, allowrepeat,
                                 dempeningFactor)
        
        # apply new results
        population[object,] = mutResult$newGenome
        evalVals[object] = NA;
        mutationCount = mutationCount + 1;
      }
      verbose(paste(mutationCount, "mutations applied\n"));
    }

    ##########
    # Add Headless Chickens
    verbose("Adding New Chromosomes ... ");

    if (newPerGen > 0) {
      for (i in (popSize+1+1):totalPopulation) { # don't mutate the parent
        population[i,] = ga.new.chromosome(genomeLen, genomeMin, genomeMax, allowrepeat)
      }
    }
    
    ##########
    # Evaluation
    verbose("Calucating evaluation values... ");
    
    to.eval.Ids = which(is.na(evalVals))
    evalVals[to.eval.Ids] = unlist(plapply(to.eval.Ids, 
                                           function(i, population, evalFunc) evalFunc(population[i, ]),
                                           population, evalFunc))
    
    # check for invalid items
    if ((!all(is.numeric(evalVals))) |
          any(is.na(evalVals)) |
          any(is.nan(evalVals))) {
      stop("Invalid cost function return value (NA or NaN).")
    }
    
    ##########
    # Selection
    # here it means only copying the best to the parent object
    verbose("  sorting results...\n");
    bestInd =  which.min(evalVals)
    #if (length(bestInd) > 1) bestInd = bestInd[1]
    
    parent  = population[bestInd,]
    parentEval = evalVals[bestInd]

    ##########
    # extract statistics about generation
    bestEvals[iter] = min(evalVals);
    meanEvals[iter] = mean(evalVals);
    verbose(" done.\n");
    
    collect.results <- function() {
      settings = list(genomeMin=genomeMin, genomeMax=genomeMax,
                      popSize=popSize, newPerGen = newPerGen, totalPopulation = popSize + newPerGen,
                      iterations=iterations, suggestions=suggestion,
                      mutationChance=mutationChance)
      
      pop.info = list(population=population, evaluations=evalVals, best=bestEvals, mean=meanEvals, currentIteration=iter)
      
      best = list(genome=population[bestInd,], cost = evalVals[bestInd]);
      
      ret = list(settings = settings, population = pop.info, best = best)
      
      class(ret) = "EvolutionStrategy.int";
      return (ret)
    }
    
    if (!is.null(monitorFunc)) {
      verbose("Sending current state to the monitor()...\n");
      monitorFunc(collect.results());
    }
    
    ##########
    # check termination conditions
    if (iter == iterations) {
      verbose("End of generations iteration reached.\n");
      break
    }
    
    if (!is.na(terminationCost)) {
      if (parentEval <= terminationCost) {
        verbose("Cost better than termination cost reached.\n");
        break
      }
    }
  }    

  # report the results
  return(collect.results());
}
