genZ2Z <-
function (genZ) 
{
    ans <- NULL
    nloc <- ncol(genZ)/3
    for (i in 1:nrow(genZ)) {
        t <- 1
        for (l in 1:nloc) {
            t <- as.vector(genZ[i, (3 * l - 2):(3 * l)]) %x% 
                t
        }
        ans <- append(ans, t)
    }
    ans <- matrix(ans, nrow = nrow(genZ), byrow = TRUE)
    rownames(ans) <- rownames(genZ)
    colnames(ans) <- genNames(nloc)
    return(ans)
}
