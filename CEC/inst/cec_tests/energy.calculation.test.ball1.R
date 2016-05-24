testname <- "Energy calculation (ball1)"
setup <- function()
{
    B <- as.matrix(read.table(system.file("cec_tests", "ball1.data", package="CEC")))
    C <- as.matrix(read.table(system.file("cec_tests", "centers1.data", package="CEC")))
}


test.type.covariance <- function()
{
    given.cov = matrix(c(2,1,1,3), 2,2)  
    
    expected.energy <- 2.766927173
    
    CE <- cec(B, centers=1, type="cov", param = given.cov, iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.fixedr <- function()
{
    r <- 1.5
    
    expected.energy <- 2.410818718
    
    CE <- cec(B, centers=1, type="fix", param = 1.5, iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.spherical <- function()
{
    expected.energy <- 1.456430201
    
    CE <- cec(B, centers=1, type="sp", iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}


test.type.diagonal <- function()
{
    cov <- cov.mle(B)
    
    expected.energy <- 1.45637452
    
    CE <- cec(B, centers=1, type="diag", iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.eigenvalues <- function()
{
    evals <- c(0.1, 0.22)
    
    expected.energy <- 1.734310397
    
    CE <- cec(B, centers=1, type="eigen", param=evals, iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.all <- function()
{
    expected.energy <- 1.455903678
    
    CE <- cec(B, centers=1, type="all", iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

####################################################################################################################

test.type.spherical.cluster.removing <- function()
{
    expected.energy <- 1.456430201
    
    CE <- cec(B, C, type="sp", iter.max=20)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations], msg="Energy")
    
}


