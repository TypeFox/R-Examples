preparelinearGPmap <-
function (gmap, reference = "F2", freqmat = NULL, S_full = TRUE) 
{
    if (sum(is.na(gmap)) > 0) {
        stop("Error: vector of genotypic values contains missing values")
    }
    ans <- list()
    ans$gmap <- as.matrix(gmap)
    ans$nloc <- log(length(ans$gmap), 3)
    if (S_full) {
        S <- freqmat2Sgenofreq(ans$nloc, reference, freqmat)
        ans$smat <- S$smat
        ans$sinv <- S$sinv
        ans$genofreq <- S$genofreq
        ans$genofreqloc <- S$genofreqloc
    }
    else {
        ans$smat <- NULL
        ans$sinv <- NULL
        ans$genofreq <- 1
        ans$genofreqloc <- NULL
        for (l in 1:ans$nloc) {
            loc <- freqmat2Sgenofreqloc(reference = reference, 
                l, freqmat)
            ans$smat <- c(ans$smat, list(loc$smat))
            ans$sinv <- c(ans$sinv, list(loc$sinv))
            ans$genofreq <- kronecker(loc$genofreq, ans$genofreq)
            ans$genofreqloc <- rbind(ans$genofreqloc, loc$genofreq)
        }
    }
    return(ans)
}
