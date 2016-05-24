processDataMutSurv <- function (data.mut, data.surv, colTime=2, colStatus=3) {
    #############################################################################################
    # Find matched cases in the two data sets and perform survival analysis on gene pairs 
    # in matched cases.
    # 
    # Input:
    # data.mut: col 1: gene; other cols: one case per col and each an integer vector with three
    #           values (0 wildtype, 1 amplification, -1 deletion).
    # data.surv: matrix of data frame containing case ID, survival time and survival status.  
    #            Cases do not have to match data.mut.
    # Output:
    # same as computeSurvivalPValue.
    #############################################################################################
    
    # find matching cases in data.mut and data.surv
    data.mut.t <- t (data.mut[,-1])
    colnames (data.mut.t) <- data.mut[,1]
    
    cases.mut <- colnames (data.mut)[-1]
    # replace - in case IDs with .
    # so that case IDs have same format in both data sets
    cases.surv.fix <- sapply(1:length(data.surv[,1]), function(i)gsub("-", ".", data.surv[i,1]))
    
    # find indices of matched cases
    accounted.mut <- match (intersect (cases.mut, cases.surv.fix), cases.mut)
    accounted.surv <- match (intersect (cases.mut, cases.surv.fix), cases.surv.fix)
    
    # extract data for matched cases
    data.mut.fix <- data.mut.t[accounted.mut,]
    data.surv[,1] <- cases.surv.fix
    data.surv.fix <- data.surv[accounted.surv,]
    data.surv.fix[,colStatus] <- as.numeric (data.surv.fix[,colStatus]=="DECEASED")
    
    return (list (data.mut=data.mut.fix, data.surv=data.surv.fix))
}

