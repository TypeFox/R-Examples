RUV1 <-
function (Y, eta, ctl) 
{
    if (is.null(eta)) 
        return(Y)
    if (length(eta) == 1) 
        if (eta == 1) 
            eta = matrix(1, 1, ncol(Y))
    Yc = Y[, ctl, drop = FALSE]
    etac = eta[, ctl, drop = FALSE]
    return(Y - Yc %*% t(etac) %*% solve(etac %*% t(etac)) %*% 
        eta)
}
