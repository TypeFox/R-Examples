#Construct a fuzzy system from a vector of real values. (chromosome)
fugeR.constructFuzzyFromChromosome <-
function(   chromosome=c(),
            nbRule,
            nbMaxVarInPerRule,
            nbVarOut,
            nbInputSet,
            nbOutputSet        ) {

    #check chromosome length
    chromoExpectedLength <- ((nbRule-1)*((nbMaxVarInPerRule*(2+(nbInputSet))) +
                            (nbVarOut*2))) + (nbVarOut*(1+nbOutputSet))
    if (length(chromosome) != chromoExpectedLength) {
        print("Bad chromosome size, expected :")
        print(((nbRule-1)*((nbMaxVarInPerRule*(2+(nbInputSet))) +
              (nbVarOut*2))) + (nbVarOut*(1+nbOutputSet)))
        print("Actual size is :")
        print(length(chromosome))
        stop("program stopped --bad chromosome size--")
    }

    #Start decoding the chromosome
    #Remove the default rule part
    #TO DO
    chromosomeRules <- chromosome[1:(chromoExpectedLength-(nbVarOut*(1+nbOutputSet)))]
    chromosomeDefaultRule <- chromosome[(length(chromosomeRules)+1):chromoExpectedLength]
    
    #reconstruct the fuzzy system
    matSystem <- matrix(chromosomeRules, nrow=nbRule-1, byrow=T)

    #Decompose input and output
    inputMatrix <- matSystem[ , 1:(ncol(matSystem) - (nbVarOut*2)), drop=F]
    outputMatrix <- matSystem[ , (ncol(inputMatrix) + 1):ncol(matSystem), drop=F]
    
    #Find input Var ID
    colToSelect <- seq(1, ncol(inputMatrix), nbInputSet+2)
    varInIds <- apply(inputMatrix[ ,colToSelect, drop=F], c(1,2), as.integer)
    
    #Find input MF ID
    colToSelect <- seq(2, ncol(inputMatrix), nbInputSet+2)
    varInMfIds <- apply(inputMatrix[ ,colToSelect, drop=F], c(1,2), as.integer)

    #Find input MF
    seqMf <- rep(1:nbInputSet, nbMaxVarInPerRule)
    toAdd <- rep(seq(2,ncol(inputMatrix),2+nbInputSet), each=nbInputSet)
    colToSelect <- seqMf + toAdd
    varInMfs <- inputMatrix[ ,colToSelect, drop=F]

    #Find output Var ID
    colToSelect <- seq(1, ncol(outputMatrix), 2)
    varOutIds <- apply(outputMatrix[ ,colToSelect, drop=F], c(1,2), as.integer)

    #Find output MF ID
    colToSelect <- seq(2, ncol(outputMatrix), 2)
    varOutMfIds <- apply(outputMatrix[ ,colToSelect, drop=F], c(1,2), as.integer)

    #Find output var MF
    seqMf <- rep(1:nbOutputSet, nbVarOut)
    toAdd <- rep(seq(1,length(chromosomeDefaultRule),1+nbOutputSet), each=nbOutputSet)
    colToSelect <- seqMf + toAdd
    varOutMfs <- matrix(chromosomeDefaultRule[colToSelect], ncol=nbOutputSet, byrow=T)

    #Find default rule
    colToSelect <- seq(1, length(chromosomeDefaultRule), 1+nbOutputSet)
    defautMfIds <- sapply(chromosomeDefaultRule[colToSelect], as.integer)


    fuzzy = list(   
      #type="fuzzySystem",
      inputVarIds  = varInIds,
      inputMfIds   = varInMfIds,
      inputMfs     = varInMfs,
      outputVarIds = varOutIds,
      outputMfIds  = varOutMfIds,
      outputMfs    = varOutMfs,
      defautMfIds  = defautMfIds,
      nbOut        = nbVarOut,
      nbMaxIn      = nbMaxVarInPerRule,
      nbInMf       = nbInputSet,
      nbOutMf      = nbOutputSet,
      nbRule       = nbRule
    )
    
    class(fuzzy) = "fuzzySystem"
    return(fuzzy)
} 
