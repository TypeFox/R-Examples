global.1m.analysis <- function(XC, XT, A, alpha = 0.05, n = NULL) {

    if (missing(XC)) stop("Missing 'XC' argument.")
    if (missing(XT)) stop("Missing 'XT' argument.")
    if (!is.numeric(n) && !is.null(n)) stop("The 'n' argument should be numeric.")
    if ((alpha < 0) || (alpha > 1)) stop("The 'alpha' argument shuld be between 0 and 1.")
    if (!is.numeric(alpha)) stop("The 'alpha' argument should be numeric.")
    XC <- as.matrix(XC)
    XT <- as.matrix(XT)
    if (is.null(n)) n <- nrow(XC) else {if (n != nrow(XC)) stop("Wrong value of 'n'.")}
    
    # Line below Equation (10) page 385 in our 2014 paper
    Y <- rbind(XC, XT)

    m <- ncol(Y)
    
    one.vec <- rep(1.0, 2 * n)
    g.vec <- c(rep(0.0, n), rep(1.0, n))

    if (missing(A)) {
        C <- matrix(c(0, 1), nrow = 1)
        Gamma <- cbind(one.vec, g.vec)
    } else {
        if (is.matrix(A)) {
            C <- matrix(c(0, 1, rep(0, ncol(A))), nrow = 1)
            Gamma <- cbind(one.vec, g.vec, A)
        } else {
            C <- matrix(c(0, 1, 0), nrow = 1)
            Gamma <- cbind(one.vec, g.vec, A)
        }
    }

    # MLE estimation
    tol <- 10E-40
    tmp <- t(Gamma) %*% Gamma
    if (isTRUE(all.equal(det(tmp), 0))) {
	stop ("Maximum Likehood Estimation is not possible with this dataset.")
    } else {
        GammatGamma.inv <- solve(t(Gamma) %*% Gamma, tol = tol)
	Bhat <- GammatGamma.inv %*% (t(Gamma) %*% Y)
    } 


    # Computation of unbiased covariance matrix estimate
    if (missing(A)) {
        coef <- (1 / (2 * n - 2))
    } else {
        coef <- (1 / (2 * n - ncol(Gamma)))
    }
    tmp <- Y - Gamma %*% Bhat
    Sigma <- coef * t(tmp) %*% tmp


    # Computation of the p-value
    Im <- diag(1.0, m)
    A1 <- kronecker(C, Im)
    A2 <- kronecker(GammatGamma.inv, Sigma)
    A3 <- kronecker(t(C), Im)
    # Equation (11) on page 386 in our 2014 paper:
    W <- A1 %*% A2 %*% A3
    # Equation of Z_n^2 or T_n^2 on page 386 in our 2014 paper:
    tmp <- C %*% Bhat
    Zn2 <- tmp %*% solve(W) %*% t(tmp)
    pvalueZn2 <- 1 - pchisq(Zn2, m)
    
    result <- list("Pvalue" = pvalueZn2)
    
    return(result)

}
