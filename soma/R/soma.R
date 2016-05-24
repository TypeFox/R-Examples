soma <- function (costFunction, bounds, options = list(), strategy = "all2one", ...)
{
    if (!(all(c("min","max") %in% names(bounds))))
        report(OL$Error, "Bounds list must contain \"min\" and \"max\" vector elements")
    if (length(bounds$min) != length(bounds$max))
        report(OL$Error, "Bounds are not of equal length")
    if (strategy != "all2one")
        report(OL$Error, "Only the \"all2one\" strategy is currently supported")
    
    # Use defaults for unspecified options
    defaultOptions <- list(pathLength=3, stepLength=0.11, perturbationChance=0.1, minAbsoluteSep=0, minRelativeSep=1e-3, nMigrations=20, populationSize=10)
    defaultsNeeded <- setdiff(names(defaultOptions), names(options))
    spuriousOptions <- setdiff(names(options), names(defaultOptions))
    options[defaultsNeeded] <- defaultOptions[defaultsNeeded]
    if (length(spuriousOptions) > 0)
        report(OL$Warning, "The following options were specified but are not used: ", paste(spuriousOptions,collapse=", "))
    
    nParams <- length(bounds$min)
    nParamsTotal <- nParams * options$populationSize
    steps <- seq(0, options$pathLength, options$stepLength)
    nSteps <- length(steps)
    steps <- rep(steps, each=nParamsTotal)
    
    # Create the population
    population <- matrix(runif(nParamsTotal), nrow=nParams, ncol=options$populationSize)
    population <- population * (bounds$max-bounds$min) + bounds$min
    
    # Calculate initial costs
    costFunctionValues <- apply(population, 2, costFunction, ...)
    
    migrationCount <- 0
    leaderCostHistory <- numeric(0)
    
    report(OL$Info, "Starting SOMA optimisation")
    
    repeat
    {
        # Find the current leader
        leaderIndex <- which.min(costFunctionValues)
        leaderValue <- costFunctionValues[leaderIndex]
        separationOfExtremes <- max(costFunctionValues) - leaderValue
        sumOfExtremes <- max(costFunctionValues) + leaderValue
        
        # Check termination criteria
        if (migrationCount == options$nMigrations)
        {
            report(OL$Info, "Migration limit (", options$nMigrations, ") reached - stopping")
            break
        }
        if (separationOfExtremes < options$minAbsoluteSep)
        {
            report(OL$Info, "Absolute cost separation (", signif(separationOfExtremes,3), ") is below threshold (", signif(options$minAbsoluteSep,3), ") - stopping")
            break
        }
        # isTRUE() needed here in case extremes are infinite: Inf/Inf => NaN
        if (isTRUE(separationOfExtremes/sumOfExtremes < options$minRelativeSep))
        {
            report(OL$Info, "Relative cost separation (", signif(separationOfExtremes/sumOfExtremes,3), ") is below threshold (", signif(options$minRelativeSep,3), ") - stopping")
            break
        }
        
        leaderCostHistory <- c(leaderCostHistory, leaderValue)
        
        # Find the migration direction for each individual
        directionsFromLeader <- apply(population, 2, "-", population[,leaderIndex])
        
        # Establish which parameters will be changed
        toPerturb <- runif(nParamsTotal) < options$perturbationChance
        
        # Second line here has a minus because directions are away from leader
        populationSteps <- array(rep(population,nSteps), dim=c(nParams,options$populationSize,nSteps))
        populationSteps <- populationSteps - steps * rep(directionsFromLeader * toPerturb, nSteps)
        
        # Replace out-of-bounds parameters with random valid values
        outOfBounds <- which(populationSteps < bounds$min | populationSteps > bounds$max)
        randomSteps <- array(runif(nParamsTotal*nSteps), dim=c(nParams,options$populationSize,nSteps))
        randomSteps <- randomSteps * (bounds$max-bounds$min) + bounds$min
        populationSteps[outOfBounds] <- randomSteps[outOfBounds]
        
        # Values over potential locations
        costFunctionValues <- apply(populationSteps, 2:3, costFunction, ...)
        individualBestLocs <- apply(costFunctionValues, 1, which.min)
        
        # Migrate each individual to its best new location, and update costs
        indexingMatrix <- cbind(seq_len(options$populationSize), individualBestLocs)
        population <- t(apply(populationSteps, 1, "[", indexingMatrix))
        costFunctionValues <- costFunctionValues[indexingMatrix]
        
        migrationCount <- migrationCount + 1
        if (migrationCount %% 10 == 0)
            report(OL$Verbose, "Completed ", migrationCount, " migrations")
    }
    
    report(OL$Info, "Leader is #", leaderIndex, ", with cost ", signif(costFunctionValues[leaderIndex],3))
    
    returnValue <- list(leader=leaderIndex, population=population, cost=costFunctionValues, history=leaderCostHistory, migrations=migrationCount)
    class(returnValue) <- "soma"
    
    return (returnValue)
}

plot.soma <- function (x, y = NULL, ...)
{
    plot(seq_along(x$history), x$history, xlab="Migration number", ylab="Leader cost value", type="b", ...)
}
