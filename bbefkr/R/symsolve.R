symsolve <-
function(Asym, Bmat)
{
    n <- ncol(Asym)
    if(max(abs(Asym - t(Asym)))/max(abs(Asym)) > 1e-10)
		stop("Argument not symmetric.")
    Lmat <- chol(Asym, T)
    if(attr(Lmat, "rank") < n)
        stop("Argument singular.")
    Lmatinv <- solve(Lmat[, order(attr(Lmat, "pivot"))])
    Xmat <- Lmatinv %*% t(Lmatinv) %*% Bmat
    return(Xmat)
}
