####Return the fuzzified matrix value [row, col] [nbRule, nbMaxVar]
# Return a list of matrix each matrix corespond to one case in the dataset
# 
# fuzzySystemIn : represensation of a fuzzy system
# dataset       : the data to fuzzify

fugeR.fuzzify <-
function(fuzzySystem, dataset) {

    fugeR.nbMaxVarInPerRule <- fuzzySystem$nbMaxIn
    fugeR.nbInputSet <- fuzzySystem$nbInMf

    #Find which var are used by the system in the dataset
    nbCase <- nrow(dataset)
    nbVar <- ncol(dataset)
    nbRule <- fuzzySystem$nbRule
    lstVarUsed <- unique(fuzzySystem$inputVarIds[fuzzySystem$inputVarIds %in% 1:nbVar])
    
    #Rule to compute
    #Minus 1, because the default rule is computing differently
    nbRuleToCompute <- 1:(nbRule-1)
    defaultMinValue <- rep(2.0,nbCase)
    lstRule <- vector("list", (nbRule-1))

    #Compute Rules activation
    for(i in nbRuleToCompute) {
        minValue <- defaultMinValue
        for(j in 1:fugeR.nbMaxVarInPerRule) {
            idVar <- fuzzySystem$inputVarIds[i,j]
            #if var not used... we skip to next var
            if(!(idVar %in% lstVarUsed)) {
            	next
            }
            #find the mf functions of the variable
            lstMf <-  fuzzySystem$minIn[idVar] + (
		                  (fuzzySystem$inputMfs[i,
		                  ((j*fugeR.nbInputSet)-(fugeR.nbInputSet-1)):(j*fugeR.nbInputSet)]) *
	                    fuzzySystem$intervalIn[idVar])

            #***********************************************************#
            #Fuzzify the variable and apply min operator (AND)
            minValue <- cfuzzifyVar(fuzzySystem$inputMfIds[i,j], lstMf, dataset[,idVar], minValue)
            #***********************************************************#
        }
        lstRule[[i]] <- minValue
    }
    #return the evaluation for each rule
    return(lstRule)
}
