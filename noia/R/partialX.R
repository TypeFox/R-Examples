partialX <-
function (genZ, reference = "F2", effect) 
{
    loci <- which(strsplit(effect, "")[[1]] != noia::effectsNames[1])
    if (length(loci) == 0) {
        loci <- 1
        ans <- matrix(rep(1, nrow(genZ)), dimnames = list(NULL, 
            "."))
    }
    else {
        colZ <- sort(c(3 * loci - 2, 3 * loci - 1, 3 * loci))
        zs <- genZ2ZS(reference = reference, genZ = genZ[, colZ])
        ans <- (zs$zmat %*% zs$smat)
    }
    ans <- as.matrix(ans)
    decomp <- sapply(strsplit(colnames(ans), ""), c)
    if (is.vector(decomp)) {
        decomp <- t(decomp)
    }
    colnames(ans) <- apply(decomp, 2, effNames, loci, nchar(effect))
    return(ans)
}
