permcoefs.plsRglm <- function(dataset, ind, nt, modele, family = NULL, maxcoefvalues,ifbootfail,verbose){
    tempcoefs <- try(PLS_glm_wvc(dataY = dataset[ind, 1], dataX = dataset[, 
        -1], nt = nt, modele = modele, family=family, keepstd.coeffs = TRUE, verbose=verbose)$std.coeffs, silent=TRUE)
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
