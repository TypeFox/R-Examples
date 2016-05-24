upset.incidence <-
function(z, Q=NULL, ...) {
    n <- ifelse(is.logical(Q), sum(Q), length(Q))
    #as.vector(rowSums(z[,Q])>0) # problema con Q scalare
    res <- colSums(matrix(z[Q,], nrow=n))>0
    names(res) <- rownames(z)
    return(res)
}
