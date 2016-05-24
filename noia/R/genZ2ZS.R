genZ2ZS <-
function (genZ, reference = "F2", max.level = NULL, max.dom = NULL, 
    threshold = 0) 
{
    "strrev" <- function(ss) {
        sapply(lapply(strsplit(ss, character(0)), rev), paste, 
            collapse = "")
    }
    ans <- list()
    N <- nrow(genZ)
    ans$smat <- 1
    ans$zmat <- as.matrix(rep(1, N))
    nloc <- ncol(genZ)/3
    for (l in 1:nloc) {
        eff <- colnames(ans$smat)
        geno <- colnames(ans$zmat)
        ans$zmat <- t(apply(cbind(genZ[, (3 * l - 2):(3 * l)], 
            ans$zmat), 1, function(x) {
            c(x[1] * x[4:length(x)], x[2] * x[4:length(x)], x[3] * 
                x[4:length(x)])
        }))
        ans$smat <- kronecker(Sloc(reference = reference, l, 
            genZ), ans$smat)
        if (is.null(eff)) {
            colnames(ans$smat) <- noia::effectsNames[1:3]
        }
        else {
            colnames(ans$smat) <- strrev(kronecker(noia::effectsNames[1:3], 
                strrev(eff), "paste", sep = ""))
        }
        if (is.null(geno)) {
            colnames(ans$zmat) <- noia::genotypesNames
        }
        else {
            colnames(ans$zmat) <- strrev(kronecker(noia::genotypesNames, 
                strrev(geno), "paste", sep = ""))
        }
        rownames(ans$smat) <- colnames(ans$zmat)
        useful.effects <- effectsSelect(nloc = nloc, max.level = max.level, 
            max.dom = max.dom, effects = colnames(ans$smat))
        useful.genotypes <- colnames(ans$zmat)
        ans$smat <- ans$smat[useful.genotypes, useful.effects]
        ans$zmat <- ans$zmat[, useful.genotypes]
    }
    rownames(ans$smat) <- colnames(ans$zmat)
    rownames(ans$zmat) <- rownames(genZ)
    return(ans)
}
