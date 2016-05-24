dropAttr <-
function (mat) {
    d <- dim(mat)
    mat <- as.vector(mat)
    dim(mat) <- d
    mat
}
