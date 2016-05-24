BLUP <- function(model,RE=NULL) {
    ## Note on ordered - each of REmodelmatrices, Vi and z have been ordered
    if(length(RE)==0) RE <- names(model$sampinfo$REmodelmatrices)
    if(length(RE)==0) return(cat("There are no Random Effects in this model.\n"))
    if(!all(RE %in% model$args$CorModels)) return(cat("These Random Effects are not in this model.\n"))

    ## conditional mean of each RE assuming each starts life as a
    ## zero-mean RE with diagonal covariance matrix. We insert the
    ## estimate of the parsill.
    Wymf <- model$estimates$Vi %*% (model$sampinfo$z - model$sampinfo$X %*% model$estimates$betahat) ## W times (observed - fitted)
    Us <- list()
    theta <- model$estimates$theta
    type <- attributes(theta)$type
    terms <- attributes(theta)$terms

    for(ii in 1:length(RE)) {
        ind <- match(RE[ii],terms)
        ind2 <- match(RE[ii],names(model$sampinfo$REmodelmatrices))
        Us[[ii]] <- as.vector(theta[ii] * crossprod(model$sampinfo$REmodelmatrices[[ind2]], Wymf))
        names(Us[[ii]]) <- colnames(model$sampinfo$REmodelmatrices[[ind2]])
    }
    names(Us) <- RE

    ## Conditional Covariance matrices
    ks <- unlist(lapply(model$sampinfo$REmodelmatrices,ncol))
    COV <- diag(rep(theta[match(RE,terms)],ks))
    TT <- NULL
    for(ii in 1:length(RE)) TT <- cbind(TT,theta[match(RE[ii],terms)] * model$sampinfo$REmodelmatrices[[RE[ii]]])
    COV <- COV - t(TT) %*% model$estimates$Vi %*% TT

    indices <- list()
    start <- c(1,1:cumsum(ks))[1:length(ks)]
    end <- cumsum(ks)
    for(ii in 1:length(start)) indices[[ii]] <- start[ii]:end[ii]

    names(start) <- RE
    names(end) <- RE

    means <- unlist(Us)
    rownames(COV) <- colnames(COV) <- names(means)
    vars <- diag(COV)

    allInd <- NULL
    for(ii in 1:length(RE)) {
        allInd <- c(allInd,start[ii]:end[ii])
    }

    output <- list("Mean"=means[allInd],
                   "Variance"=vars[allInd],
                   "Covariance"=COV[allInd,allInd])

    return(output)

}
