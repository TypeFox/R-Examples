freqmat2Sgenofreq <-
function (nloc, reference = "F2", freqmat = NULL, sinv = TRUE) 
{
    "strrev" <- function(ss) {
        sapply(lapply(strsplit(ss, character(0)), rev), paste, 
            collapse = "")
    }
    ans <- list()
    ans$smat <- 1
    if (sinv) {
        ans$sinv <- 1
    }
    ans$genofreq <- 1
    ans$genofreqloc <- NULL
    for (l in 1:nloc) {
        eff <- colnames(ans$smat)
        geno <- rownames(ans$smat)
        loc <- freqmat2Sgenofreqloc(reference = reference, l, 
            freqmat)
        ans$smat <- kronecker(loc$smat, ans$smat)
        if (sinv) {
            ans$sinv <- kronecker(loc$sinv, ans$sinv)
        }
        if (is.null(eff)) {
            colnames(ans$smat) <- noia::effectsNames[1:3]
        }
        else {
            colnames(ans$smat) <- strrev(kronecker(noia::effectsNames[1:3], 
                strrev(eff), "paste", sep = ""))
        }
        if (is.null(geno)) {
            rownames(ans$smat) <- noia::genotypesNames
        }
        else {
            rownames(ans$smat) <- strrev(kronecker(noia::genotypesNames, 
                strrev(geno), "paste", sep = ""))
        }
        ans$genofreq <- kronecker(loc$genofreq, ans$genofreq)
        ans$genofreqloc <- rbind(ans$genofreqloc, loc$genofreq)
    }
    colnames(ans$sinv) <- rownames(ans$smat)
    rownames(ans$sinv) <- colnames(ans$smat)
    colnames(ans$genofreqloc) <- noia::genotypesNames
    colnames(ans$genofreq) <- colnames(sinv)
    return(ans)
}
