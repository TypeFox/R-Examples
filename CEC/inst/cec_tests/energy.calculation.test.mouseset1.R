testname <- "Energy calculation (mouseset1)"
setup <- function()
{
    B <- as.matrix(read.table(system.file("cec_tests", "mouse1.data", package="CEC")))
}

test.type.covariance <- function()
{
    given.cov = matrix(c(2,1,1,3), 2,2)
    
    expected.energy <- 3.540174056 
    
    CE <- cec(B, centers=1, type="cov", param = given.cov, iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.fixedr <- function()
{
    r <- 1.5
    
    expected.energy <- 3.416637007
    
    CE <- cec(B, centers=1, type="fix", param = 1.5, iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.spherical <- function()
{
    expected.energy <- 3.403158062  
    
    CE <- cec(B, centers=1, type="sp", iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}


test.type.diagonal <- function()
{
    expected.energy <- 3.396500695
    
    CE <- cec(B, centers=1, type="diag", iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}

test.type.all <- function()
{
    expected.energy <- 3.396472329
    
    CE <- cec(B, centers=1, type="all", iter.max=0)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[1], msg="Energy")
}
