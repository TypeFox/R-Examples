# 15) Compute the AV using the IF.
# Input
# IFiv=Sample Influence function of IV estimator, ((p+1)xn) matrix.
AVif <- function(IFriv) {
    n <- ncol(IFriv)
    A <- apply(IFriv, 2, tcrossprod)
    dim(A) <- c(nrow(IFriv), nrow(IFriv), n)
    
    A.AV <- apply(A, 1:2, mean)/n
    dimnames(A.AV) <- list(rownames(IFriv), rownames(IFriv))

    A.AV
}
