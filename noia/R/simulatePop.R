simulatePop <-
function (gmap, N = 100, sigmaE = 1, type = "F2", freqmat = NULL) 
{
    nloc <- log(length(gmap))/log(3)
    ans <- NULL
    n <- genNames(nloc)
    if (length(n) != length(gmap)) {
        stop("Genotype-phenotype map: not the expected size.")
    }
    if (is.null(names(gmap)) || !all((sort(n) == sort(names(gmap))))) {
        warning("The labels of the genotype-phenotype map do not fit with the expected genotypes. Relabelling.")
        names(gmap) <- n
    }
    for (i in 1:N) {
        g <- drawGenotype(nloc, type, freqmat)
        indiv <- c(sigmaE * rnorm(1) + gmap[g], unlist(strsplit(g, 
            character(0))))
        ans <- rbind(ans, as.numeric(indiv))
    }
    coln <- "phen"
    for (l in 1:nloc) {
        coln <- c(coln, paste("Loc", l, sep = ""))
    }
    colnames(ans) <- coln
    return(as.data.frame(ans))
}
