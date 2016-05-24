calculateMdist <-
function(dataset, treatment, matchVars){
    cat("Calculating Mahalanobis distances...\n")
    ## calculate the inverse covariance matrix
    icv <- solve(cov(dataset[, matchVars]))
    ## get the names of the treated
    trtnms <- row.names(dataset)[as.logical(dataset[[treatment]])]
    ## and the names of the control units
    ctlnms <- row.names(dataset)[!as.logical(dataset[[treatment]])]
    ## calculate the mahalanobis distances if not specified
    mdist <- outer(trtnms, ctlnms, FUN = mahalDist, inv.cov = icv, data = dataset)
    dimnames(mdist) <- list(trtnms, ctlnms)
    return(mdist)
}
