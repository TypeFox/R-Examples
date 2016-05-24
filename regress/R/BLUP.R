BLUP <- function(model,RE=NULL) {
    ## model - output from regress

    ## RE - vector of names of random effects, should match names of
    ## model$sigma

    ## Returns the conditional mean, variance and variance covariance
    ## matrices of random effects given the data. Default is all
    ## random effects ignoring those associated with the identity
    ## matrix.

    if(length(RE)==0) RE <- setdiff(model$Vnames,"In")
    if(any(is.na(match(RE,model$Vnames)))) stop(paste("RE should be a subset of",model$Vnames))
    if(class(model)!="regress") stop("model should be of class regress")

    ## conditional expected values - all of them
    Wy <- model$W %*% (model$model[[1]] - model$fitted)
    Us <- list()
    for(ii in 1:length(model$Z)) {
        Us[[ii]] <- as.vector(model$sigma[ii] * t(model$Z[[ii]]) %*% Wy)
        names(Us[[ii]]) <- colnames(model$Z[[ii]])
    }
    names(Us) <- names(model$sigma)

    ## Conditional Covariance matrices - ignore In for now
    ks <- unlist(lapply(model$Z,ncol))
    COV <- diag(rep(model$sigma,ks))
    TT <- NULL
    for(ii in (1:length(model$Z))) TT <- cbind(TT,model$sigma[ii] * model$Z[[ii]])
    COV <- COV - t(TT) %*% model$W %*% TT

    indices <- list()
    start <- c(1,1+cumsum(ks))[1:length(ks)]
    end <- cumsum(ks)
    for(ii in 1:length(start)) indices[[ii]] <- start[ii]:end[ii]

    names(start) <- names(model$sigma)
    names(end) <- names(model$sigma)

    means <- unlist(Us)
    rownames(COV) <- colnames(COV) <- names(means)
    vars <- diag(COV)

    allInd <- NULL
    for(ii in 1:length(RE)) {
        ind <- match(RE[ii],model$Vnames)
        ##output[[ii]] <- list("Mean"=means[start[ind]:end[ind]],
        ##                     "Variance"=vars[start[ind]:end[ind]],
        ##                     "Covariance"=COV[start[ind]:end[ind],start[ind]:end[ind]])
        allInd <- c(allInd,start[ind]:end[ind])
    }
    ##names(output) <- RE
    output <- list("Mean"=means[allInd],
                   "Variance"=vars[allInd],
                   "Covariance"=COV[allInd,allInd])
    return(output)
}
