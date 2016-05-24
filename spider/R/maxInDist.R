maxInDist <- 
function(distobj, sppVector = NULL, propZero = FALSE, rmNA = FALSE){
    dat <- as.matrix(distobj)
    if(length(sppVector) > 0) dimnames(dat)[[1]] <- sppVector
    conSpecDists <- list()
    for (i in 1:length(dimnames(dat)[[1]])) {
        conSpec <- dimnames(dat)[[1]] == dimnames(dat)[[1]][i]
        conSpecDists[[i]] <- max(dat[conSpec, i], na.rm = rmNA)
    }
    if (propZero) 
        output <- length(which(unlist(conSpecDists) == 0))/length(unlist(conSpecDists))
    else output <- unlist(conSpecDists)
    output
}
