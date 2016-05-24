computeSurvivalPValueForGenePairSet.output <- function (file.out, gene.pairs, data.mut, data.surv, colTime=2, colStatus=3, type.gene1=(-1), type.gene2=(-1), groups=c("All", "Two"), PRINT=FALSE, PRINT.INDEX=FALSE) {
    #############################################################################################
    # Find matched cases in the two data sets and perform survival analysis on a given set of 
    # gene pairs in matched cases.  Write to file directly.
    # 
    # Input:
    # data.mut: col 1: gene; other cols: one case per col and each an integer vector with three
    #           values (0 wildtype, 1 amplification, -1 deletion).
    # data.surv: matrix of data frame containing case ID, survival time and survival status.  
    #            Cases do not have to match data.mut.
    # Output:
    # same as computeSurvivalPValue.
    #############################################################################################
    groups = match.arg (groups)
    
    # process mutation and survival data
    data <- processDataMutSurv (data.mut=data.mut, data.surv=data.surv, colTime=colTime, colStatus=colStatus)
    data.mut.fix <- data$data.mut
    data.surv.fix <- data$data.surv
    
    # remove NAs in survival data
    which.na <- which (is.na (data.surv.fix[,colTime]) | is.na (data.surv.fix[,colStatus]))
    data.mut <- data.mut.fix[-which.na,]
    data.surv <- data.surv.fix[-which.na,]
    
    gene.pairs.new <- data.frame (gene1=gene.pairs[,1], gene1mut=type.gene1, gene2=gene.pairs[,2], gene2mut=type.gene2)    

    # run survival analysis
    cat ("computing p values...\n")
    computeSurvivalPValueGenePairAll.output (file.out, gene.pairs.new, data.mut, data.surv, colTime=colTime, colStatus=colStatus, groups=groups, PRINT=PRINT, PRINT.INDEX=PRINT.INDEX)

}

