procrustus <- function (x) {
    sx <- svd (x)
    return (tcrossprod (sx$u, sx$v))
}