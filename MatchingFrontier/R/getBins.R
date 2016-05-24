getBins <-
function(dataset, treatment, match.on, breaks){
    gs <- getL1Strata(treatment, dataset, match.on, breaks=breaks)
    strata <- gs$strata
    mycut <- gs$mycut
    names(strata) <- dataset[,which(colnames(dataset) == treatment)]
    unique.strata <- unique(strata)
    strataholder <- list()
    for(i in 1:length(unique.strata)){
        strataholder[[i]] <- which(strata==unique.strata[i]) 
    }
    return(strataholder)
}
