matSums <-
function (lis) {
    out <- array(data = 0, dim = dim(lis[[1]]))
    for (i in seq(along = lis))
        out <- out + lis[[i]]
    out
}
