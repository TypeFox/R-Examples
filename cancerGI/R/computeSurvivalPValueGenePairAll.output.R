computeSurvivalPValueGenePairAll.output <- function (file.out, gene.pairs, data.mut, data.surv, colTime=2, colStatus=3, groups=c("All", "Two"), PRINT=FALSE, PRINT.INDEX=FALSE) {
    groups = match.arg (groups)
    
    if (groups=="All") {
        pqvalues <- data.frame (nNoMut=0, nGene1Mut=0, nGene2Mut=0, nDoubleMut=0, obsNoMut=0, obsGene1Mut=0, obsGene2Mut=0, obsDoubleMut=0, expNoMut=0, expGene1Mut=0, expGene2Mut=0, expDoubleMut=0, medianNoMut=0, medianGene1Mut=0, medianGene2Mut=0, medianDoubleMut=0, pValue=NA)        
    }
    else {
        pqvalues <- data.frame (nSingleMut=0, nDoubleMut=0, obsSingleMut=0, obsDoubleMut=0, expSingleMut=0, expDoubleMut=0, medianSingleMut=0, medianDoubleMut=0, pValue=NA)
    }
    
    gene1.id <- match (gene.pairs[,1], colnames (data.mut))
    gene2.id <- match (gene.pairs[,3], colnames (data.mut))
    
    pqvalues[1,] <- computeSurvivalPValueOneGenePair (data.mut[,c(gene1.id[1], gene2.id[1])], data.surv, colTime=colTime, colStatus=colStatus, type.gene1=gene.pairs[1,2], type.gene2=gene.pairs[1,4], groups=groups, PRINT=PRINT)[-1]
    
    genes <- data.frame (gene.pairs, gene1ID=gene1.id, gene2ID=gene2.id)
    write (c(colnames(genes)[1:4], colnames(pqvalues)), file.out, sep="\t", ncolumns=4+ncol(pqvalues))
    
    #apply (gene.pairs[-1,], 1, testfn.output, file.out=file.out, append=TRUE)
    apply (genes, 1, computeSurvivalPValueOneGenePair.output, file.out=file.out, data.mut=data.mut, data.surv=data.surv, colTime=colTime, colStatus=colStatus, groups=groups, PRINT=PRINT)
}

