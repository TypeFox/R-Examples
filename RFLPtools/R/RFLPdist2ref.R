###############################################################################
## Computation of distances for RFLP data
###############################################################################

## x: data.frame with RFLP data
## ref: data.frame with RFLP reference data
## distfun: function to compute distance (cf. ?dist)
RFLPdist2ref <- function(x, ref, distfun = dist, nrBands, LOD = 0){
    stopifnot(is.data.frame(x))
    stopifnot(is.data.frame(ref))
    stopifnot(is.function(distfun))
    
    if(missing(nrBands))
        stop("Number of Bands 'nrBands' is missing.")
    if(nrBands <= 0)
        stop("'nrBands' has to be a positive interger!")
    if(LOD > 0){
        x <- x[x$MW >= LOD,]
        ref <- ref[ref$MW >= LOD,]
    }
        
    x1 <- split(x, x$Sample)
    ref1 <- split(ref, ref$Sample)
    nrbands <- sort(unique(sapply(x1, nrow)))
    ref.nrbands <- sort(unique(sapply(ref1, nrow)))
    
    if(!(nrBands %in% nrbands))
        stop("There is no sample with specified number of bands 'nrBands'.")
        
    if(!(nrBands %in% ref.nrbands))
        stop("There is no reference sample with specified number of bands 'nrBands'.")
     
    x1.bands <- sapply(x1, nrow)
    ref1.bands <- sapply(ref1, nrow)

    temp <- do.call("rbind", x1[x1.bands == nrBands])
    ref.temp <- do.call("rbind", ref1[ref1.bands == nrBands])

    temp1 <- split(temp[,"MW"], factor(temp[,"Sample"]))
    grp <- factor(paste(ref.temp[,"Taxonname"], " (", ref.temp[,"Accession"], ")", sep = ""))
    ref.temp1 <- split(ref.temp[,"MW"], grp)

    res <- as.matrix(distfun(do.call("rbind", c(ref.temp1, temp1))))
    res[(length(ref.temp1)+1):nrow(res), 1:length(ref.temp1), drop = FALSE]
}
