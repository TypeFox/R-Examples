permcoefs.plsR <- function(dataset,ind,nt,modele,maxcoefvalues,ifbootfail,verbose){
tempcoefs <- try(PLS_lm_wvc(dataY =dataset[ind,1], dataX=dataset[,-1], nt=nt, keepstd.coeffs=TRUE, verbose=verbose)$std.coeffs,silent=TRUE)
    Cond <- FALSE
    try(Cond<-is.numeric(tempcoefs)&all(abs(tempcoefs)<maxcoefvalues),silent=TRUE)
    if (Cond) {
        return(tempcoefs)
    }
    else {
        return(ifbootfail)
    }
}
