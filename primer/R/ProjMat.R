`ProjMat` <-
function (stagedat, fruitdat, seeddat) 
{
    fecs <- tapply(fruitdat$Y2004, fruitdat$Stage, mean)/2
    seed.freqs <- table(seeddat[, 1])
    seedfates <- seed.freqs/length(seeddat[, 1])
    seedfates
    mat1 <- matrix(0, nrow = 5, ncol = 5)
    for (i in 2:5) {
        for (j in 2:5) mat1[i, j] <- {
            x <- subset(stagedat, stagedat$Y2003 == j)
            jT <- nrow(x)
            iT <- sum(x$Y2004 == i)
            iT/jT
        }
    }
    mat1[1, 1] <- seedfates[2]
    mat1[2, 1] <- seedfates[3]
    mat1[1, 4] <- fecs[1]
    mat1[1, 5] <- fecs[2]
    return(mat1)
}
