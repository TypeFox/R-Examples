getL1Strata <-
function(treatment, dataset, match.on, breaks){
    # Remove dropped covs
    dataset <- dataset[,match.on]
    ## stuff borrowed from cem.main to add user defined breaks
    vnames <- colnames(dataset)
    nv <- dim(dataset)[2]
    mycut <- vector(nv, mode="list")
    names(mycut) <- vnames
    for (i in 1:nv) {	
        tmp <- reduceVar(dataset[[i]], breaks[[vnames[i]]])
        dataset[[i]] <- tmp$x
        mycut[[vnames[i]]] <- tmp$breaks
    }
    # Calculate strata
    strata <- stratify(dataset)
    return(list(strata=strata, mycut=mycut))
}
