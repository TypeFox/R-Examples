#Return the fuzzified matrix value [row, col] [1, nbOut]
#     outputFuzzy <- list(
#                             varIdMat = matrix(as.integer(fuzzySystem$outputVarIds), nbRule),
#                             mfIdMat  = matrix(as.integer(fuzzySystem$outputMfIds), nbRule),
#                             mfValMat = fuzzySystem$outputMfs
#                        )
fugeR.defuzzify <- function(fuzzySystem, fuzzifiedValue) {
    
    nbVarOut    <- fuzzySystem$nbOut
    nbRule      <- fuzzySystem$nbRule

    #Find var and mf
    lstInRule <- vector( 'list', nbVarOut )
    lstMfId   <- vector( 'list', nbVarOut )
    lstMf     <- vector( 'list', nbVarOut )
    
    #Find membership functions
    sapply(1:nbVarOut, function(x) {
      lstInRule[[x]] <<- fuzzySystem$outputVarIds[ ,x] %in% 1:nbVarOut
      lstMfId[[x]]   <<- fuzzySystem$outputMfIds[ ,x]
      lstMf[[x]]     <<- sort(fuzzySystem$minOut[x] + (fuzzySystem$outputMfs[x,] *
                         fuzzySystem$intervalOut[x]))
    } )
    
    prediction <- vector( 'list', nbVarOut )
    
    #This defuzzify output var
    for( i in 1:nbVarOut ) {
      prediction[[i]] <- cdefuzzify(  length(fuzzifiedValue[[1]]),
                                      nbRule,
                                      lstInRule[[i]],
                                      lstMf[[i]],
                                      lstMfId[[i]],
                                      fuzzySystem$defautMfIds[i],
                                      fuzzifiedValue)
    }
    
    return(prediction)
} 
