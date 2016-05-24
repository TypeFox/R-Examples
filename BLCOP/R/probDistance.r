probDistance <- function(result, numSimulations = BLCOPOptions("numSimulations")) 
{
    monteCarloSample <- rmnorm(numSimulations, result@priorMean, result@priorCovar)
    mean(abs(dmnorm(monteCarloSample, result@priorMean, result@priorCovar,log=TRUE) - dmnorm(monteCarloSample, result@posteriorMean, 
        result@posteriorCovar, log = TRUE)))
}

setGeneric("probDistance")

probDistance.COPResult <- function(result, numSimulations = BLCOPOptions("numSimulations") )
{
    show("Not implemented yet...")
}

setMethod("probDistance", signature(result = "COPResult"), probDistance.COPResult)