## compute stress per point

spp <- function(dhat, confdiss, wgths) 
{
    resmat <- as.matrix(wgths)*as.matrix(dhat - confdiss)^2    #point stress
    diag(resmat) <- NA
    spp <- colMeans(resmat, na.rm = TRUE)
    spp <- spp/sum(spp)*100
    names(spp) <- colnames(resmat) <- rownames(resmat) <- attr(dhat, "Labels") 
    return(list(spp = spp, resmat = resmat))
}

