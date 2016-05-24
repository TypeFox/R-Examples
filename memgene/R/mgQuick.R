mgQuick <-
function(genD, coords, longlat=FALSE, truncation=NULL, transformation=NULL,
                    forwardPerm=100, forwardAlpha=0.05, finalPerm=NULL,
                    doPlot=NULL, verbose=TRUE) {
    
    if (longlat) {
        if (verbose) cat("Finding geodesic distances among sampling locations\n")
        if (!require(geosphere)) {
                    stop("memgene: geosphere package must be installed to use longlat=TRUE", call.=FALSE)
        }
        eucDist <- as.matrix(distm(coords, fun=distMeeus))
        eucDist[is.na(eucDist)] <- 0
    }
    else {
        if (verbose) cat("Finding Euclidean distances among sampling locations\n")
        eucDist <- as.matrix(dist(coords))
    }
    
    if (verbose) cat("Extracting Moran's eigenvectors\n")
    mem <- tryCatch(mgMEM(eucDist, truncation=truncation, transformation=transformation),
                    error=function(e) return(list(error=e,
                                            P=NA,
                                            RsqAdj=0,
                                            F=NA,
                                            memgene=NA,
                                            memSelected=NA,
                                            whichSelectedPos=NA,
                                            whichSelectedNeg=NA,
                                            resid=NA,
                                            pred=NA,
                                            sdev=NA,
                                            mem=NA)))
    if (!is.null(mem$error)) return(mem)
    pos <- mem$vectorsMEM[, mem$valuesMEM > 0]
    neg <- mem$vectorsMEM[, mem$valuesMEM < 0]

    if (verbose) cat("Forward selections of positive Moran's eigenvectors\n")
    selectedPos <- tryCatch(na.omit(mgForward(genD, pos, perm=forwardPerm,
                    alpha=forwardAlpha)$selectedMEM),
                      error=function(e) return(list(error=e,
                                            P=NA,
                                            RsqAdj=0,
                                            F=NA,
                                            memgene=NA,
                                            memSelected=NA,
                                            whichSelectedPos=NA,
                                            whichSelectedNeg=NA,
                                            resid=NA,
                                            pred=NA,
                                            sdev=NA,
                                            mem=mem)))
    if (class(selectedPos) == "list") return(selectedPos)
    if (verbose) cat("----Selected:", ifelse(length(selectedPos)==0, "None", paste(sort(selectedPos), collapse=", ")), "\n")
    
    if (verbose) cat("Forward selections of negative Moran's eigenvectors\n")
    selectedNeg <- tryCatch(na.omit(mgForward(genD, neg, perm=forwardPerm,
                        alpha=forwardAlpha)$selectedMEM),
                            error=function(e) return(list(error=e,
                                            P=NA,
                                            RsqAdj=0,
                                            F=NA,
                                            memgene=NA,
                                            memSelected=NA,
                                            whichSelectedPos=NA,
                                            whichSelectedNeg=NA,                                      
                                            resid=NA,
                                            pred=NA,
                                            sdev=NA,
                                            mem=mem)))
    if (class(selectedNeg) == "list") return(selectedNeg)
    if (verbose) cat("----Selected:", ifelse(length(selectedNeg)==0, "None", paste(sort(selectedNeg), collapse=", ")), "\n")  
    memSelected <- cbind(pos[, selectedPos], neg[, selectedNeg])
    if (ncol(memSelected)==0) {
        if (verbose) cat("There are no significant Moran's eigenvectors\n")
        return(list(error=simpleError("No significant Moran's eigenvectors"),
            P=NA,
            RsqAdj=0,
            F=NA,
            memgene=NA,
            memSelected=NA,
            whichSelectedPos=NA,
            whichSelectedNeg=NA,              
            resid=NA,
            pred=NA,
            sdev=NA,
            mem=mem))
    }
    else {
        if (verbose) cat("Finding memgene variables using selected Moran's eigenvectors\n")    
        if (verbose && !is.null(finalPerm)) cat("Running permutation test\n")
        final <- tryCatch(mgRDA(genD, memSelected, perm=finalPerm),
                          error=function(e) return(list(error=e,
                                            P=NA,
                                            RsqAdj=0,
                                            F=NA,
                                            memgene=NA,
                                            memSelected=memSelected,
                                            whichSelectedPos=selectedPos,
                                            whichSelectedNeg=selectedNeg,
                                            resid=NA,
                                            pred=NA,
                                            sdev=NA,
                                            mem=mem)))
        if (!is.null(final$error)) return(final)                            
        
        if (!is.null(doPlot)) {
            if (verbose) cat("Mapping first", doPlot, "memgene variable(s)\n")
            mgMap(coords, final$memgene[, 1:doPlot])
        }
        
        result <-list(P=final$P,
            RsqAdj=final$RsqAdj,
            F=final$F,
            memgene=final$memgene,
            memSelected=memSelected,
            whichSelectedPos=selectedPos,
            whichSelectedNeg=selectedNeg,            
            resid=final$resid,
            pred=final$pred,
            sdev=final$sdev,
            mem=mem)
        
        class(result) <- "mgQuick"
        return(result)
    }
}
