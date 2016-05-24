mvsf <-
function (m) 
{
    if (!is.matrix(m)) 
        stop("m[] is not a matrix with number of columns (sample size) between 5 and 5000")
    n <- ncol(m)
    if (n < 5 || n > 5000) 
        stop("sample size must be between 5 and 5000")
    R <- m - apply(m, 1, mean)
    M.1 <- solve(R %*% t(R), tol = 1e-18)
    Rmax <- diag(t(R) %*% M.1 %*% R)
    C <- M.1 %*% R[, which.max(Rmax)]
    Z <- t(C) %*% m
    a<-sf.test(Z)$statistic
    b<-sf.test(Z)$p.value
    c<-sf.test(Z)$data.name
    results <- list(statistic = a, p.value = b, method = "Generalized Shapiro-Francia test for Multivariate Normality", data.name = c)
    class(results) = "htest"
    return(results)
}

