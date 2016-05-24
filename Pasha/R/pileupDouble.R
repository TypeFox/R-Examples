

pileupDouble <- function(start, fragLength, dir,readLength, weight=1) 
{
    
    # The c code assume that readlength, fraglength and dir are specified for all reads
    nbReads <- length(start)
    
    # Recycle values for all reads in case a unique value has been provided (can be optimised, see martin and simon method)
    if(length(fragLength)==1) fragLength <- rep(fragLength, nbReads)
    if(length(dir)==1) dir <- rep(dir, nbReads)
    if(length(readLength)==1) readLength <- rep(readLength, nbReads)
    if(length(weight)==1) weight <- rep(weight, nbReads)
    
    
    dir <- factor(dir, levels=c("+", "-"))
    
    # drop non-standard values for strand direction
    nonStandard <- is.na(dir)
    
    if(any(nonStandard))
    {
        warning("Some undefined strand (dir) values were ignored...")
        start <- start[!nonStandard]
        fragLength <- fragLength[!nonStandard]
        dir <- dir[!nonStandard]
        readLength <- readLength[!nonStandard]
        weight <- weight[!nonStandard]
    }
    
    stopifnot((length(start)==nbReads), (length(fragLength)==nbReads), (length(dir)==nbReads), (length(readLength)==nbReads), (length(weight)==nbReads))
    
    # Compute the max coordinate possible (theoretical) in order to create in mem the resulting vector of double
    maxCoord <- max(start)+max(c(max(fragLength), max(readLength)))
    
    functionReturn <- .C("C_pileupDouble", as.integer(start), as.integer(fragLength), as.integer(dir), as.integer(readLength), as.double(weight), as.integer(nbReads), as.integer(maxCoord), res=double(maxCoord))
    
    return(functionReturn$res)
}



