prepareRegression <-
function (phen, gen = NULL, genZ = NULL, reference = "noia", 
    max.level = NULL, max.dom = NULL, fast = FALSE) 
{
    ans <- list()
    ans$phen <- phen
    ans$gen <- gen
    ans$genZ <- genZ
    ans$phen <- as.matrix(ans$phen)
    ans$gen <- as.matrix(ans$gen)
    if (!is.null(ans$gen)) {
        if (nrow(ans$phen) != nrow(ans$gen)) {
            stop("Error: not the same number of genotypes and phenotypes.")
        }
        ans$gen <- as.matrix(ans$gen[!is.na(ans$phen), ])
    }
    if (!is.null(ans$genZ)) {
        if (nrow(ans$phen) != nrow(ans$genZ)) {
            stop("Error: not the same number of genotypes and phenotypes.")
        }
        ans$genZ <- ans$genZ[!is.na(ans$phen), ]
        checkgenZ(ans$genZ)
    }
    else {
        if (is.null(ans$gen)) {
            stop("No genotype provided.")
        }
        else {
            ans$genZ <- gen2genZ(ans$gen)
        }
    }
    ans$phen <- ans$phen[!is.na(ans$phen)]
    ans$nloc <- ncol(ans$genZ)/3
    ans$reference <- reference
    SandZ <- list()
    ans$x <- NULL
    if (fast) {
        ans$x <- genZ2X(genZ = ans$genZ, reference = reference, 
            max.level = max.level, max.dom = max.dom)
    }
    else {
        SandZ <- genZ2ZS(genZ = ans$genZ, reference = reference, 
            max.level = max.level, max.dom = max.dom)
        ans$zmat <- SandZ$zmat
        ans$smat <- SandZ$smat
        ans$x <- ans$zmat %*% ans$smat
    }
    return(ans)
}
