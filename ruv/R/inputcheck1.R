inputcheck1 <-
function (Y, X, Z, ctl) 
{
    m = nrow(Y)
    n = ncol(Y)
    if (m > n) 
        warning("m is greater than n!")
    if (sum(is.na(Y)) > 0) 
        warning("Y contains missing values.  This is not supported.")
    if (sum(Y == Inf) + sum(Y == -Inf) > 0) 
        warning("Y contains infinities.  This is not supported.")
    XZ = cbind(X, Z)
    d = svd(t(XZ) %*% XZ, nu = 0, nv = 0)$d
    if (d[1]/d[length(d)] > 10^10) 
        warning("There appears to be linear dependence between the columns of X and Z.")
    if ((length(ctl) != n) | (!is.logical(ctl))) 
        warning("ctl should be a logical vector of length n.")
    if (sum(ctl) == 0) 
        warning("No genes are defined as control genes.  This is not supported.")
    return(NULL)
}
