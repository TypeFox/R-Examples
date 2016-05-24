 #Evalutate a fuzzy system during evolution
fugeR.evaluate <-
function(string=c()) {
    #Parameters
    fitness = NA
    fugeR.nbRule <- fugeRglobal[["fugeR.nbRule"]]
    fugeR.nbMaxVarInPerRule <- fugeRglobal[["fugeR.nbMaxVarInPerRule"]]
    fugeR.nbVarOut <- fugeRglobal[["fugeR.nbVarOut"]]
    fugeR.nbInputSet <- fugeRglobal[["fugeR.nbInputSet"]]
    fugeR.nbOutputSet <- fugeRglobal[["fugeR.nbOutputSet"]]
    
    #Data
    fugeR.ACTUAL_VALUES <- fugeRglobal[["fugeR.ACTUAL_VALUES"]]
    fugeR.TO_PREDICT <- fugeRglobal[["fugeR.TO_PREDICT"]]
    
    #Construct the fuzzy from the chromosome
    fuzzy <- fugeR.constructFuzzyFromChromosome(  string,
                                                  fugeR.nbRule,
                                                  fugeR.nbMaxVarInPerRule,
                                                  fugeR.nbVarOut,
                                                  fugeR.nbInputSet,
                                                  fugeR.nbOutputSet )
    
    fuzzy$intervalIn  = fugeRglobal$fugeR.intervalIn
    fuzzy$minIn       = fugeRglobal$fugeR.minIn
    fuzzy$intervalOut = fugeRglobal$fugeR.intervalOut
    fuzzy$minOut      = fugeRglobal$fugeR.minOut
    
    #Make the prediction on the training dataset (the one provided)
    predictedValues <- fugeR.predict(fuzzy, fugeR.ACTUAL_VALUES)
    
    #Calculate the fitness of the fuzzySystem and return it
    fitness <- fugeR.calcFitness(predictedValues, fugeR.TO_PREDICT, fuzzy)
    return(fitness)
}
