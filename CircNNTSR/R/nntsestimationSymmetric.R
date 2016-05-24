nntsestimationSymmetric<-
function (M = 0, data, maxit = 500) 
{
    if (M == 0) 
        return("Uniform distribution: no parameters to estimate")
    A <- rbind(diag(1, M + 1, M + 1), cbind(matrix(0, nrow = 1, 
        ncol = M), diag(-1, 1, 1)), c(rep(-1, M), rep(0, 1)))
    cvector <- rbind(matrix(0, nrow = M + 1, ncol = 1), matrix(-2 * 
        pi + 1e-11, nrow = 1, ncol = 1), -1/(2 * pi))
    IP <- nntsrandominitialSymmetric(M)
    r <- constrOptim(theta = IP[1:(M + 1)], f = nntsloglikSymmetric, 
        grad = NULL, ui = A, ci = cvector, method = "Nelder-Mead", 
        data = data, M = M, control = list(maxit = maxit))
    AIC <- (2 * r$value) + 2 * (M + 1)
    BIC <- (2 * r$value) + (M + 1) * log(length(data))
    res <- list(coef = r$par, loglik = -r$value, AIC = AIC, BIC = BIC, 
        convergence = r$convergence)
    return(res)
}
