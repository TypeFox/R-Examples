###############################################################################
## Computation of distances for RFLP data
###############################################################################

## x: data.frame with RFLP data
## distfun: function to compute distance (cf. ?dist)
RFLPdist <- function(x, distfun = dist, nrBands, LOD = 0){
    stopifnot(is.data.frame(x))
    stopifnot(is.function(distfun))
    if(!missing(nrBands)){
        if(nrBands <= 0)
            stop("'nrBands' has to be a positive interger!")
    }
    if(LOD < 0)
        stop("'LOD' has to be non-negative!")

    if(LOD > 0){
        x <- x[x$MW >= LOD,]
    }

    x1 <- split(x, x$Sample)
    nrbands <- sort(unique(sapply(x1, nrow)))
    x1.bands <- sapply(x1, nrow)

    if(missing(nrBands)){
        res <- vector("list", length(nrbands))
        index <- 0
        for(i in nrbands){
            index <- index + 1
            temp <- do.call("rbind", x1[x1.bands == i])
            temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
            res[[index]] <-  distfun(do.call("rbind", temp1))
        }
        names(res) <- nrbands
    }else{
        if(!(nrBands %in% nrbands))
            stop("No samples with given number of bands!")
        temp <- do.call("rbind", x1[x1.bands == nrBands])
        temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
        res <-  distfun(do.call("rbind", temp1))
    }
    res
}
