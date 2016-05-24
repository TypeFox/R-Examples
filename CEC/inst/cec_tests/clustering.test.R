testname <- "Clustering"
setup <- function()
{
    M <- as.matrix(read.table(system.file("cec_tests", "mouse1.data", package="CEC")))
    C <- as.matrix(read.table(system.file("cec_tests", "centers1.data", package="CEC")))
    expected <- dget(system.file("cec_tests", "cec1.dp", package="CEC"))  
}
test.clustering.mouse1 <- function()
{  
    CEC <- cec(M, C)
    CEC:::checkNumericVectorEquals(expected$cluster, CEC$cluster, msg="Clustering vector")
    CEC:::checkNumericVectorEquals(expected$cost, CEC$cost, msg="Energy")
    CEC:::checkNumericMatrixEquals(expected$centers, CEC$centers, msg="Centers")
    CEC:::checkNumericMatrixEquals(M, CEC$data, msg="Data")
    CEC:::checkNumericVectorEquals(expected$probability, CEC$probability, msg="Probability")
    CEC:::checkNumericVectorEquals(expected$nclusters, CEC$nclusters, msg="Number of clusters")
    CEC:::checkNumericVectorEquals(expected$final.cost.function, CEC$final.cost.function, msg="Final cost function")
    CEC:::checkNumericVectorEquals(expected$final.nclusters, CEC$final.nclusters, msg="Final number of clusters")
    CEC:::checkNumericVectorEquals(expected$iterations, CEC$iterations, msg="Iterations")
    
}
