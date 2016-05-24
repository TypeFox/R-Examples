GeneticAlg.int <- function(genomeLen, codonMin, codonMax,
                           genomeMin=rep.int(codonMin, genomeLen), genomeMax=rep.int(codonMax, genomeLen),
                           suggestions=NULL,
                           popSize=50, 
                           iterations=100, terminationCost=NA,
                           mutationChance=1/(genomeLen+1),
                           elitism=floor(popSize/10),
                           geneCrossoverPoints=NULL,
                           monitorFunc=NULL, 
                           evalFunc,
                           allowrepeat = TRUE,
                           showSettings=FALSE, verbose=FALSE,
                           plapply = lapply) {
  # Optimizes an Integer chromosome using a genetic algorithm.
  #
  # popSize          = the population size
  # iterations       = number of generations
  # terminationCost  = The cost (error) that if reached, the GA should termiante
  # mutationChance   = chance that a var in the string gets mutated
  # geneCrossoverPoints = An array determining the genes to be swapped in crossover
  #
  # Partially based on "R Based Genetic Algorithm (genalg package)""
  # http://cran.r-project.org/web/packages/genalg/
  
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
  if (popSize < 5) {
    stop("The population size must be at least 5.");
  }
  if (iterations < 1) {
    stop("The number of iterations must be at least 1.");
  }
  if (!(elitism < popSize)) {
    stop("The population size must be greater than the elitism.");
  }
  if (elitism < 0) {
    stop("elitism must be at least 0.");
  }
  if ((mutationChance < 0) | (mutationChance  > 1)) {
    stop("mutationChance must be between 0 and 1.");
  }
  if (!is.null(geneCrossoverPoints)) {
    if (!is.numeric(geneCrossoverPoints) | length(geneCrossoverPoints) != genomeLen) {
      stop("Invalid geneCrossoverPoints.");
    }
  }
  
  if (showSettings) {
    verbose("The start conditions:\n");
    result = list(genomeMin=genomeMin, genomeMax=genomeMax, suggestions=suggestions,
                  popSize=popSize, iterations=iterations,
                  elitism=elitism, mutationChance=mutationChance);
    class(result) = "rbga";
    
    cat(summary(result));
  } else {
    verbose("Not showing GA settings...\n");
  }
  
  ##########
  # Creation
  population = matrix(nrow=popSize, ncol=genomeLen);
  
  if (!is.null(suggestions)) {
    verbose("Adding suggestions to first population...\n");
    suggestionCount = nrow(suggestions)
    population[1:suggestionCount,] = population[i,]
    verbose("Filling others with random values in the given domains...\n");
  } else {
    verbose("Starting with random values in the given domains...\n");
    suggestionCount = 0
  }
  
  for (i in (suggestionCount+1):popSize) {
    population[i,] = ga.new.chromosome(genomeLen, genomeMin, genomeMax, allowrepeat)
  }
  
  ############################################################################
  # do iterations
  bestEvals = rep(NA, iterations);
  meanEvals = rep(NA, iterations);
  evalVals = rep(NA, popSize);
  for (iter in 1:iterations) {
    verbose(paste("Starting iteration", iter, "\n"));
    
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
    
    # extract statistics about generation
    bestEvals[iter] = min(evalVals);
    meanEvals[iter] = mean(evalVals);
    bestInd = which.min(evalVals)
    verbose(" done.\n");

    collect.results <- function() {
      settings = list(genomeMin=genomeMin, genomeMax=genomeMax,
                      popSize=popSize, elitism=elitism, geneCrossoverPoints = geneCrossoverPoints,
                      iterations=iterations, suggestions=suggestions,
                      mutationChance=mutationChance)
      
      pop.info = list(population=population, evaluations=evalVals, best=bestEvals, mean=meanEvals, currentIteration=iter)
      
      best = list(genome=population[bestInd,], cost = evalVals[bestInd]);
      
      ret = list(settings = settings, population = pop.info, best = best)
      
      class(ret) = "GeneticAlg.int";
      return (ret)
    }
    
    if (!is.null(monitorFunc)) {
      verbose("Sending current state to rgba.monitor()...\n");
      # report on GA results
      monitorFunc(collect.results());
    }
    
    ##########
    # check termination conditions
    if (iter == iterations) {
      verbose("End of generations iteration reached.\n");
      break
    }
    
    if (!is.na(terminationCost)) {
      if (bestEvals[iter] <= terminationCost) {
        verbose("Cost better than termination cost reached.\n");
        break
      }
    }
    
    ##########
    # Selection
    
    verbose("Creating next generation...\n");
    newPopulation = matrix(nrow=popSize, ncol=genomeLen);
    newEvalVals = rep(NA, popSize);
    
    verbose("  sorting results...\n");
    sortedEvaluations = sort(evalVals, index=TRUE);
    sortedPopulation  = matrix(population[sortedEvaluations$ix,], ncol=genomeLen);
    
    # save the best
    if (elitism > 0) {
      verbose("  applying elitism...\n");
      newPopulation[1:elitism,] = sortedPopulation[1:elitism,];
      newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
    } # ok, save nothing
    
    ##########
    # Crossover
    # fill the rest by doing crossover
    verbose("  applying crossover...\n");
    for (child in (elitism+1):popSize) {
      # ok, pick two random parents using roulette wheel probability
      parentProb = dnorm(1:popSize, mean=0, sd=(popSize/3))
      parentIDs = sample(1:popSize, 2, prob=parentProb)
      parents = sortedPopulation[parentIDs,]
      
      # default cross-over swaps genomes from a random point
      if (is.null(geneCrossoverPoints)) {
        crossOverPoint = sample(0:genomeLen,1)
        
        if (crossOverPoint == 0) {
          newPopulation[child, ] = parents[2,]
          newEvalVals[child] = sortedEvaluations$x[parentIDs[2]]
        } else if (crossOverPoint == genomeLen) {
          newPopulation[child, ] = parents[1,]
          newEvalVals[child] = sortedEvaluations$x[parentIDs[1]]
        } else {
          newPopulation[child, ] = 
            c(parents[1, 1:crossOverPoint], 
              parents[2, (crossOverPoint+1):genomeLen])
        }
      } else {
        # swap genes based on the location of gene crossoverPoints
        p2.indices = which(geneCrossoverPoints %% 2 != 0)
        
        newPopulation[child, ] = parents[1,]
        newPopulation[child, p2.indices] = parents[2, p2.indices]
      }
    }
    
    if (!allowrepeat) {
      for (i in (elitism+1):popSize) {
        population[i,] = ga.unique.maker(population[i,], genomeMin, genomeMax)
      }
    }
    
    population = newPopulation;
    evalVals   = newEvalVals;
    
    ##########
    # Mutation
    if (mutationChance > 0) {
      verbose("  applying mutations... ");
      mutationCount = 0;
      for (object in (elitism+1):popSize) { # don't mutate the best
        
        dempeningFactor = (iterations-iter)/iterations
        
        mutResult <- ga.mutation(population[object,], mutationChance, genomeLen, 
                                 genomeMin, genomeMax, allowrepeat,
                                 dempeningFactor)
        
        population[object, ] = mutResult$newGenome;
        evalVals[object] = NA;
        mutationCount = mutationCount + 1;
      }
      verbose(paste(mutationCount, "mutations applied\n"));
    }
  }
  
  # report on GA results
  result = collect.results()
  
  return(result);
}
