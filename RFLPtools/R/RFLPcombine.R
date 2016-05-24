###############################################################################
## Combine data sets
###############################################################################

RFLPcombine <- function(...){
    x <- list(...)
    if(length(x) < 2)
        stop("You have to specify at least two data sets!")
    nams0 <- unique(x[[1]]$Sample)
    for(i in 2:length(x)){
        n0 <- length(nams0)
        ## extract sample names for dataset i
        nams1 <- unique(x[[i]]$Sample)
        ## make unique names
        nams0 <- make.unique(c(nams0, nams1))
        ## extract unique names for dataset i
        nams2 <- nams0[(n0+1):length(nams0)]
        ## names that have been changed
        nams3 <- nams1[nams1 != nams2]
        nams4 <- nams2[nams1 != nams2]
        ## replace names that have been changed by unique names
        for(j in 1:length(nams3)){
            x[[i]]$Sample[x[[i]]$Sample == nams3[j]] <- nams4[j]
        }
    }
    do.call('rbind', x)
}
