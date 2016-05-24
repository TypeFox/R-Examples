norm = function (A, p = 2)
{
    A <- as.matrix(A)
    if (min(dim(A)) == 1)
        A <- t(A)
    if (p == 1)
        return(as.matrix(max(colSums(abs(A)))))
    else if (p == 2) {
        A.sv <- La.svd(A)$d
        return(as.matrix(max(A.sv)))
    }
    else if (p > 1e+09)
        return(as.matrix(max(rowSums(abs(A)))))
    else stop("Unknown norm")
}
