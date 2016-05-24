diffDist <- function(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2){ 
    x.diff <- t(diff(t(as.matrix(x)))) 
    res <- dist(x.diff, method = method, diag = diag, upper = upper, p = p)
    attr(res, "call") <- match.call()
    res
}
