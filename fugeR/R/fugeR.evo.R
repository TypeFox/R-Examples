# It optimizes a binary chromosome using a genetic algorithm.
#
# string           = string to optimize
# popSize          = the population size
# iters            = number of generations
# mutationChance   = chance that a var in the string gets mutated
fugeR.evo <- function(stringMin=c(), stringMax=c(),
                 suggestions=NULL,
                 popSize=200, iters=100, 
                 mutationChance=NA,
                 elitism=NA,
                 monitorFunc=NULL, evalFunc=NULL,
                 showSettings=FALSE, verbose=FALSE) {
    if (is.null(evalFunc)) {
        stop("A evaluation function must be provided. See the evalFunc parameter.");
    }
    
    vars = length(stringMin);
    if (is.na(mutationChance)) {
        mutationChance = 1/(vars+1);
    }
    if (is.na(elitism)) {
        elitism = floor(popSize/5)
    }
    
    # TODO: should do a variaty of sanity checks first
#     if (verbose) cat("Testing the sanity of parameters...\n");
    if (length(stringMin) != length(stringMax)) {
        stop("The vectors stringMin and stringMax must be of equal length.");
    }
    if (popSize < 5) {
        stop("The population size must be at least 5.");
    }
    if (iters < 1) {
        stop("The number of iterations must be at least 1.");
    }
    if (!(elitism < popSize)) {
        stop("The population size must be greater than the elitism.");
    }
    
    if (showSettings) {
        if (verbose) cat("The start conditions:\n");
        result = list(stringMin=stringMin, stringMax=stringMax, suggestions=suggestions,
                      popSize=popSize, iters=iters,
                      elitism=elitism, mutationChance=mutationChance);
        class(result) = "rbga";
        
        cat(summary(result));
    } else {
#         if (verbose) cat("Not showing GA settings...\n");
    }
    
    if (vars > 0) {
        if (!is.null(suggestions)) {
            if (verbose) cat("Adding suggestions to first population...\n");
            population = matrix(nrow=popSize, ncol=vars);
            suggestionCount = dim(suggestions)[1]
            for (i in 1:suggestionCount) {
                population[i,] = suggestions[i,]
            }
            if (verbose) cat("Filling others with random values in the given domains...\n");
            for (var in 1:vars) {
                population[(suggestionCount+1):popSize,var] = stringMin[var] +
                                   runif(popSize-suggestionCount)*(stringMax[var]-stringMin[var]);
            }
        } else {
            if (verbose) cat("Starting with random values in the given domains...\n");
            # start with an random population
            population = matrix(nrow=popSize, ncol=vars);
            # fill values
            for (var in 1:vars) {
                population[,var] = stringMin[var] +
                                   runif(popSize)*(stringMax[var]-stringMin[var]);
            }
        }
        
        # do iterations
        bestEvals = rep(NA, iters);
        meanEvals = rep(NA, iters);
        evalVals = rep(NA, popSize);
        
        #Default value to set at each iteration
        initPopulation = matrix(nrow=popSize, ncol=vars);
        initEvalVals = rep(NA, popSize);
        
        for (iter in 1:iters) {
            if (verbose) cat(paste("GENERATION : ", iter, "\n"));
            
            haveToBeEvaluate <- is.na(evalVals)
            # calculate each object
            for (object in 1:popSize) {
                if (haveToBeEvaluate[object]) {
                    evalVals[object] = evalFunc(population[object,]);
                }
            }
            
            bestEvals[iter] = min(evalVals);
            meanEvals[iter] = mean(evalVals);
  
            if (!is.null(monitorFunc)) {
#                 if (verbose) cat("Sending current state to rgba.monitor()...\n");
                # report on GA settings
                result = list(type="floats chromosome",
                              stringMin=stringMin, stringMax=stringMax,
                              popSize=popSize, iter=iter, iters=iters,
                              population=population, elitism=elitism, mutationChance=mutationChance,
                              evaluations=evalVals, best=bestEvals, mean=meanEvals);
                class(result) = "rbga";
                
                monitorFunc(result);
            }
            #cat('\n', 'DEBUG1')
            if (iter < iters) { # ok, must create the next generation
                # if (verbose) cat("Creating next generation...\n");
                newPopulation = initPopulation;
                newEvalVals = initEvalVals;
                
                #if (verbose) cat("  sorting results...\n");
                sortedEvaluations = sort(evalVals, index=TRUE);
                sortedPopulation  = matrix(population[sortedEvaluations$ix,], ncol=vars);
                
                # save the best
                if (elitism > 0) {
                    #if (verbose) cat("  applying elitism...\n");
                    newPopulation[1:elitism,] = sortedPopulation[1:elitism,];
                    newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
                } # ok, save nothing
                #cat('\n', 'DEBUG2')
                
                #Random dad and mom
                dadIDs = sample(1:popSize,popSize-elitism)
                momIDs = sample(1:popSize,popSize-elitism)
                crossOverPoints = sample(1:(vars-1),popSize-elitism,replace=TRUE)
                coupleIndex <- 1
                # fill the rest by doing crossover
                for (child in (elitism+1):popSize) {
                  #make couple have a baby
                  newPopulation[child, ] = 
                    c(sortedPopulation[dadIDs[coupleIndex],][1:crossOverPoints[coupleIndex]], 
                      sortedPopulation[momIDs[coupleIndex],][(crossOverPoints[coupleIndex]+1):vars])
                  
                  coupleIndex <- coupleIndex + 1
                }
                
                population = newPopulation;
                evalVals   = newEvalVals;
                haveToMutate <- runif((popSize-elitism)*vars) < mutationChance
                
                # do mutation
                indexMutation <- 1
                if (mutationChance > 0) {
                    for (object in (elitism+1):popSize) { # don't mutate the best
                        for (var in 1:vars) {
                            #if (runif(1) < mutationChance) { # ok, do mutation
                          if (haveToMutate[indexMutation]) {
                                # OPTION 1
                                # mutate to something random
                                mutation = stringMin[var] +
                                   runif(1)*(stringMax[var]-stringMin[var]);
                                
                                # apply mutation, and delete known evalutation value
                                population[object,var] = mutation;
                                evalVals[object] = NA;
                            }
                          indexMutation <- indexMutation + 1
                        }
                    }
                }
            }
        }
    }

    # report on GA settings
    result = list(type="floats chromosome", 
                  stringMin=stringMin, stringMax=stringMax,
                  popSize=popSize, iters=iters, suggestions=suggestions,
                  population=population, elitism=elitism, mutationChance=mutationChance,
                  evaluations=evalVals, best=bestEvals, mean=meanEvals);
    class(result) = "rbga";

    return(result);
}

 
