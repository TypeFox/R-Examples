testMutationalPatternPair <- function (a, mtx, PRINT=FALSE) {
    return (testMutationalPatternBinom (mtx[,a[1]], mtx[,a[2]], PRINT=PRINT))
}

