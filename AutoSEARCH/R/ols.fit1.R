ols.fit1 <-
function(y, x, tol = 1e-07 , LAPACK = FALSE)
{
out <- list()
qx <- qr(x, tol, LAPACK = LAPACK)
out <- c(out, qx)
out$coefficients <- as.vector(solve.qr(qx, y, tol = tol))
return(out)
}
