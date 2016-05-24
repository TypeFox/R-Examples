downset.incidence <-
function(z, Q=NULL, ...) {
    n <- ifelse(is.logical(Q), sum(Q), length(Q))
    res <- rowSums(matrix(z[,Q], ncol=n))>0
    names(res) <- colnames(z)
    return(res)
}
