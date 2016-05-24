xy <-
function(n, p.all, p.true, R2, beta0=5, yname=paste('y', p.true, sep=''), xname=paste('x', p.true, p.all, sep='')) {
	X = matrix(rnorm(n * (p.all-1)), n)
    ey = if (p.true==1) rep(beta0, n) else beta0 + as.matrix(X[ , 1:(p.true-1)]) %*% rep(1, p.true-1)
    stdev = sqrt(crossprod(ey - beta0) * (1/R2 - 1) / n)
    err = rnorm(n, sd=stdev)
    y = ey + err
    assign(xname, X, env=.GlobalEnv)
    assign(yname, y, env=.GlobalEnv)
    list(X=X, y=y)
}

