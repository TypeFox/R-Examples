cdf.qp <-
function (expectreg, x = NA, qout = NA, extrap = FALSE) 
{
    penalty.term = 0.001
    epsilon = 1e-05
    max.iter = 20
    step.halfing = 0.5
    p = expectreg$asymmetries
    if (is.na(x)) 
        e <- expectreg$fitted[1, ]
    else if (length(which(expectreg$covariates[[1]] == x)) > 
        0) 
        e <- expectreg$fitted[which(expectreg$covariates[[1]] == 
            x)[1], ]
    else e <- expectreg$fitted[which.min(abs(expectreg$covariates[[1]] - 
        x))[1], ]
    e = sort(e)
    K <- length(e)
    if (is.null(p) == TRUE) {
        p <- seq(0 + 1/(K + 1), 1 - 1/(K + 1), length = K)
    }
    mu05 <- approx(p, y = e, xout = 0.5)$y
    delta <- rep(1/(K + 1), K)
    loop <- 0
    iter.diff <- 1
    while ((loop < max.iter) & (iter.diff > epsilon)) {
        loop <- loop + 1
        F <- cumsum(delta)
        Fs <- (kronecker(matrix(1:K), matrix(1, 1, K)) >= kronecker(matrix(1, 
            K, 1), t(matrix(1:K)))) * 1
        eg <- c(min(e) - 1e-04, e)
        eg <- eg[-1] + eg[-length(eg)]
        eg <- eg/2
        G <- cumsum(eg * delta)
        Gs <- kronecker(matrix(1, K, 1), t(matrix(eg))) * Fs
        h <- e - ((1 - p) * G + p * (mu05 - G))/((1 - p) * F + 
            p * (1 - F))
        hs <- -kronecker(matrix((1 - 2 * p)/((1 - p) * F + p * 
            (1 - F))), matrix(1, 1, K)) * Gs + kronecker(matrix((((1 - 
            p) * G + p * (mu05 - G)) * (1 - 2 * p))/(((1 - p) * 
            F + p * (1 - F))^2)), matrix(1, 1, K)) * Fs
        Ls <- t(hs) %*% h
        Lss1 <- t(hs) %*% hs
        Lss <- Lss1
        Lss <- Lss + penalty.term * diag(K)
        dvec <- -Ls
        Dmat <- Lss
        Amat <- cbind(diag(K), -matrix(1, K, 1))
        bvec <- matrix(c(-delta, sum(delta) - 1))
        xsi <- solve.QP(Dmat, dvec, Amat, bvec, meq = 0)$solution
        delta <- delta + step.halfing * xsi
        iter.diff <- max(abs(xsi))
    }
    F <- cumsum(delta)
    if (any(is.na(qout))) 
        qout = p
    if (extrap) 
        quant <- as.vector(my.approx(F, e, xout = qout, rule = 3)$y)
    else quant <- as.vector(my.approx(F, e, xout = qout, rule = 2)$y)
    result = list(x = e, density = delta, cdf = F, quantiles = quant, 
        qout = qout)
    class(result) = "expectilecdf"
    return(result)
}
