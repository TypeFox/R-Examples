ols.fit2 <-
function(y, x, tol = 1e-07 , LAPACK = FALSE)
{
tmp <- crossprod(x)
out <- qr(tmp)
out$xtxinv <- solve.qr(out, tol=tol, LAPACK=LAPACK)
out$xtx <- tmp
out$xty <- crossprod(x,y)
out$coefficients <- as.vector(out$xtxinv%*%out$xty)
return(out)
}
