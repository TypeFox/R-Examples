#' Generate a design matrix from raw RNAi data.
#' 
#' Function takes the raw RNAi data as input and generates a design matrix for regression.
#'
#' @param test Input; RNAi data.
#' @return a design matrix.
#' @export
#'

computeSurvivalPValueGenePairAll <- function (gene.pairs, data.mut, data.surv, colTime=2, colStatus=3, groups=c("All", "Two"), compare=c("Both", "Gene1", "Gene2"), PRINT=FALSE, PRINT.INDEX=FALSE) {
    groups = match.arg (groups)
    
    if (groups=="All") {
        pqvalues <- data.frame (nNoMut=0, nGene1Mut=0, nGene2Mut=0, nDoubleMut=0, obsNoMut=0, obsGene1Mut=0, obsGene2Mut=0, obsDoubleMut=0, expNoMut=0, expGene1Mut=0, expGene2Mut=0, expDoubleMut=0, medianNoMut=0, medianGene1Mut=0, medianGene2Mut=0, medianDoubleMut=0, pValue=NA)        
    }
    else {
        pqvalues <- data.frame (nSingleMut=0, nDoubleMut=0, obsSingleMut=0, obsDoubleMut=0, expSingleMut=0, expDoubleMut=0, medianSingleMut=0, medianDoubleMut=0, pValue=NA)
    }
    
    gene1.id <- match (gene.pairs[,1], colnames (data.mut))
    gene2.id <- match (gene.pairs[,3], colnames (data.mut))
    
    for (i in 1:nrow(gene.pairs)) {
        if (PRINT.INDEX) {
            cat (i," ")
        }
        #pqvalues[i,1] <- as.factor (gene.pairs[i,1])
        #pqvalues[i,2] <- as.factor (gene.pairs[i,3])
        #gene1.id <- which (colnames (data.mut)==gene.pairs[i,1])
        #gene2.id <- which (colnames (data.mut)==gene.pairs[i,3])
        if (PRINT) {
            print (gene.pairs[i,1])
            print (gene.pairs[i,3])
        }
        pqvalues[i,] <- computeSurvivalPValueOneGenePair (data.mut[,c(gene1.id[i], gene2.id[i])], data.surv, colTime=colTime, colStatus=colStatus, type.gene1=gene.pairs[i,2], type.gene2=gene.pairs[i,4], groups=groups, compare=compare, PRINT=PRINT)[-1]
    }
    
    pqvalues <- data.frame (gene.pairs, pqvalues)
    if (PRINT) {
        print (pqvalues)        
    }
    #qValue <- qvalue (pqvalues$pValue)$qvalue
    #pqvalues <- data.frame (pqvalues, qValue=qValue)
    return (pqvalues)
}

