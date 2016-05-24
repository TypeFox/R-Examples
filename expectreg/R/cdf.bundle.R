cdf.bundle <-
function (bundle, qout = NA, extrap = FALSE) 
{
    if (!inherits(bundle, "bundle") && !inherits(bundle, "restricted")) 
        stop("Function needs 'expectreg' estimated by bundle or restricted.")
    basis = bundle$design
    np = length(bundle$intercepts)
    pp <- bundle$asymmetries
    rst <- (bundle$response - (basis %*% bundle$trend.coef))/(basis %*% 
        bundle$residual.coef)
    u <- seq(1.2 * min(rst), 1.2 * max(rst), length = 100)
    m <- length(u)
    bundle$asymmetry.coef = sort(bundle$asymmetry.coef)
    A <- matrix(0, np + 1, m)
    A[np + 1, ] <- 1
    for (k in 1:np) {
        a1 <- (1 - pp[k]) * (u - bundle$asymmetry.coef[k]) * 
            (u <= bundle$asymmetry.coef[k])
        a2 <- pp[k] * (u - bundle$asymmetry.coef[k]) * (u > bundle$asymmetry.coef[k])
        A[k, ] <- a1 + a2
    }
    D <- diff(diag(m), diff = 2)
    lambda <- 100
    P <- lambda * t(D) %*% D
    v1 <- solve(t(A) %*% A + P, t(A) %*% c(rep(0, np), 1))
    q <- c(rep(0, np), 1)
    lambda2 <- 1
    D2 <- diff(diag(m), diff = 3)
    P2 <- lambda2 * t(D2) %*% D2
    z <- log(v1 - min(v1) + 0.02 * max(v1))
    for (it in 1:20) {
        g <- exp(z)
        r <- q - A %*% g
        B <- A * outer(rep(1, np + 1), as.vector(g))
        Q <- t(B) %*% B
        znew <- solve(Q + P2, t(B) %*% r + Q %*% z)
        dz <- max(abs(z - znew))
        z <- znew
        cat("iteration: ", it, ", convergence: ", dz, "\n")
        if (dz < 1e-06) 
            break
    }
    dens = g/(u[2] - u[1])
    F = cumsum(dens)
    dens = dens/max(F)
    F = F/max(F)
    if (any(is.na(qout))) 
        qout = pp
    if (extrap) 
        quant <- as.vector(my.approx(F, u, xout = qout, rule = 3)$y)
    else quant <- as.vector(my.approx(F, u, xout = qout, rule = 2)$y)
    result = list(x = u, density = dens, cdf = F, quantiles = quant, 
        qout = qout, random = rst)
    class(result) = c("expectilecdf", "bundledensity")
    result
}
