mat2list <-
function(x, MARGIN=1)
{
    n <- dim(x)[MARGIN]
    if (MARGIN==1) {
        out <- lapply(1:n, function(z) x[z,])
        names(out) <- rownames(x)
    } else {
        out <- lapply(1:n, function(z) x[,z])
        names(out) <- colnames(x)
    }
    out
}
