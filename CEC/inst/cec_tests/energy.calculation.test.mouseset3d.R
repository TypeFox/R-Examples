testname <- "Energy calculation and covariances (mouseset3d)"
setup <- function()
{
    B <- as.matrix(read.table(system.file("cec_tests", "mouse3d.data", package="CEC")))
    C <- as.matrix(read.table(system.file("cec_tests", "centers3d.data", package="CEC")))
    C4 <- as.matrix(read.table(system.file("cec_tests", "centers43d.data", package="CEC")))
}

test.type.covariance <- function()
{
    given.cov = matrix(c(0.770118878, 0.005481129, -0.005991149, 0.005481129, 0.766972716, 0.008996509, -0.005991149, 0.008996509,  0.821481768), 3, 3)
    
    expected.energy <-  4.365855156
    
    CE <- cec(B, centers=C, type="cov", param = given.cov, iter.max=20)

    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}

test.type.fixedr.mixture <- function()
{
    r <- c(0.2, 0.3, 0.4) 

    expected.energy <- 4.853461033
    
    CE <- cec(B, centers=C, type=c("fi", "fi", "fi"), param = r)
    
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}


test.type.spherical.one.cluster.removed <- function()
{  
    expected.energy <- 4.179257781
    expected.number.of.clusters <- 3
     
    CE <- cec(B, C4, type="sp")
    
    CEC:::checkNumericVectorEquals(expected.number.of.clusters, CE$final.nclusters, msg="Number of clusters")
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}

test.type.diagonal.spherical.mixture <- function()
{  
    expected.energy <- 4.177793598
    expected.number.of.clusters <- 3
    
    CE <- cec(B, C, type=c("diag", "diag", "sp"))  

    CEC:::checkNumericVectorEquals(expected.number.of.clusters, CE$final.nclusters, msg="Number of clusters")
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}

test.type.eigenvalues.all.fixedr.mixture <- function()
{ 
    evals1 <- c(0.8240634, 0.7739987, 0.7595220)
    evals2 <- c(0.7240634, 0.5739987, 0.3595220)
    r <- 1.0
    
    expected.energy <- 4.323007035
    expected.number.of.clusters <- 3
    
    CE <- cec(B, C4, type=c("all", "eigen", "fixedr", "eigen"), param=list(evals1, r, evals2))
    
    CEC:::checkNumericVectorEquals(expected.number.of.clusters, CE$final.nclusters, msg="Number of clusters")
    CEC:::checkNumericVectorEquals(expected.energy, CE$cost[CE$iterations + 1], msg="Energy")
}
