matArrays <-
function (lis, na.rm = FALSE) {
    M <- length(lis)
    dims <- dim(lis[[1]])
    mat <- unlist(lis)
    dim(mat) <- c(length(mat)/M, M)
    array(rowMeans(mat), dims)
}
