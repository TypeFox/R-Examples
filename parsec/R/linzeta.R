linzeta <-
function(lin) {
    n <- length(lin)
    res <- .C("linzeta", lin=as.integer(lin), n=as.integer(n),
        result=integer(n*n))
    res <- matrix(res$result, n, n, byrow = TRUE)
    rownames(res) <- colnames(res) <- names(lin)
    class(res) <- "incidence"
    return(res)
}
