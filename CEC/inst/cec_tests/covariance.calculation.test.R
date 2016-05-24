testname <- "Covariance calculation"
setup <- function()
{
    B <- as.matrix(read.table(system.file("cec_tests", "ball1.data", package="CEC")))
    centers <- as.matrix(read.table(system.file("cec_tests", "centers2.data", package="CEC")))
}

test.covariances.before.first.iteraion <- function()
{
    M <- matrix(c(-1, 102, 141, -1, 104, 2, -1, -1, 12, 4), 5, 2)
    cov <- cov.mle(M)
    C <- cec(M, centers=1, iter.max=0)
    checkNumericMatrixEquals(cov, C$covariances[[1]], msg="Covariances")
}

test.covariances.after.point.movements.between.clusters <- function()
{  
    cov <- cov.mle(B)
    C <- cec(B, centers=centers, type="sp")
    checkNumericMatrixEquals(cov, C$covariances[[1]], msg="Covariances")
}
