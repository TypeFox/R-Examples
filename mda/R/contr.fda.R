contr.fda <-
function (p = rep(1, d[1]), contrast.default = contr.helmert(length(p))) 
{
    d <- dim(contrast.default)
    sqp <- sqrt(p/sum(p))
    x <- cbind(1, contrast.default) * outer(sqp, rep(1, d[2] + 
        1))
    qx <- qr(x)
    J <- qx$rank
    qr.qy(qx, diag(d[1])[, seq(2, J)])/outer(sqp, rep(1, J - 
        1))
}

