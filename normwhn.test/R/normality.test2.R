normality.test2 <-
function (x) 
{
    n <- as.double(nrow(x))
    nvars <- as.double(ncol(x))
    m1 <- c(numeric(nvars))
    m2 <- c(numeric(nvars))
    m3 <- c(numeric(nvars))
    m4 <- c(numeric(nvars))
    for (p in 1:nvars) m1[p] <- mean(x[, p])
    for (p in 1:nvars) m2[p] <- (n^(-1)) * sum((x[, p] - m1[p])^2)
    for (p in 1:nvars) m3[p] <- (n^(-1)) * sum((x[, p] - m1[p])^3)
    for (p in 1:nvars) m4[p] <- (n^(-1)) * sum((x[, p] - m1[p])^4)
    sk <- c(numeric(nvars))
    k <- c(numeric(nvars))
    for (p in 1:nvars) sk[p] <- m3[p]/(m2[p]^(3/2))
    for (p in 1:nvars) k[p] <- m4[p]/m2[p]^2
    sk
    print("sk")
    print(sk)
    k
    print("k")
    print(k)
    mean.vector <- t(as.matrix(cov.wt(x)$center))
    covariance.matrix <- cov.wt(x)$cov
    diagcov <- diag(covariance.matrix)
    idiagcov <- 1/sqrt(diagcov)
    v.matrix <- matrix(data = 0, nrow = nvars, ncol = nvars)
    for (i in 1:nvars) v.matrix[i, i] <- idiagcov[i]
    correlation.matrix <- cor(x)
    lambda.matrix <- diag(eigen(correlation.matrix)$values)
    h.matrix <- eigen(correlation.matrix)$vectors
    xhat.matrix <- x - t(matrix(rep(mean.vector, n), nvars))
    xhatprime.matrix <- t(xhat.matrix)
    rprime.matrix <- h.matrix %*% solve(lambda.matrix)^(1/2) %*% 
        t(h.matrix) %*% v.matrix %*% xhatprime.matrix
    rprime <- rprime.matrix
    trprime <- t(rprime)
    v1 <- c(numeric(nvars))
    v2 <- c(numeric(nvars))
    v3 <- c(numeric(nvars))
    v4 <- c(numeric(nvars))
    for (j in 1:nvars) v1[j] <- mean(trprime[, j])
    for (j in 1:nvars) v2[j] <- (n^(-1)) * sum((trprime[, j] - 
        v1[j])^2)
    for (j in 1:nvars) v3[j] <- (n^(-1)) * sum((trprime[, j] - 
        v1[j])^3)
    for (j in 1:nvars) v4[j] <- (n^(-1)) * sum((trprime[, j] - 
        v1[j])^4)
    ac <- matrix(data = NA, nrow = n, ncol = nvars)
    for (j in 1:nvars) ac[, j] <- acf(trprime[, j], lag.max = n - 
        1, type = "covariance")$acf
    F3 <- c(numeric(nvars))
    F4 <- c(numeric(nvars))
    for (j in 1:nvars) F3[j] <- ac[1, j]^3 + 2 * sum(ac[2:n, 
        j]^3)
    for (j in 1:nvars) F4[j] <- ac[1, j]^4 + 2 * sum(ac[2:n, 
        j]^4)
    rtb1 <- c(numeric(nvars))
    b2 <- c(numeric(nvars))
    for (j in 1:nvars) rtb1[j] <- v3[j]/sqrt(F3[j])
    for (j in 1:nvars) b2[j] <- v4[j]/sqrt(F4[j])
    print("rtb1")
    print(rtb1)
    print("b2")
    print(b2)
    beta <- (3 * (n^2 + 27 * n - 70) * (n + 1) * (n + 3))/((n - 
        2) * (n + 5) * (n + 7) * (n + 9))
    w2 <- (-1) + sqrt(2 * (beta - 1))
    delta <- 1/sqrt(log(sqrt(w2)))
    f <- (w2 - 1)/2
    g <- (n + 1) * (n + 3)/(6 * (n - 2))
    h <- sqrt(f * g)
    y <- rtb1 * h
    z1 <- delta * log(y + sqrt(y^2 + 1))
    print("z1")
    print(z1)
    del <- ((n - 3) * (n + 1) * (n^2 + 15 * n - 4))
    aye <- ((n - 2) * (n + 5) * (n + 7) * (n^2 + 27 * n - 70))/(6 * 
        del)
    cee <- ((n - 7) * (n + 5) * (n + 7) * (n^2 + 2 * n - 5))/(6 * 
        del)
    alp <- aye + ((rtb1^2) * cee)
    kap <- ((n + 5) * (n + 7) * (n^3 + 37 * n^2 + 11 * n - 313))/(12 * 
        del)
    chi <- (b2 - 1 - rtb1^2) * (2 * kap)
    chi <- abs(chi)
    z2 <- (((chi/(2 * alp))^(1/3)) - 1 + (1/((9 * alp)))) * sqrt(9 * 
        alp)
    print("z2")
    print(z2)
    pvalsk <- c(numeric(nvars))
    for (p in 1:nvars) pvalsk[p] <- pnorm(z1[p], lower.tail = FALSE)
    for (p in 1:nvars) pvalsk[p] <- 2 * pvalsk[p]
    for (p in 1:nvars) if (pvalsk[p] > 1) 
        pvalsk[p] <- 2 - pvalsk[p]
    print("H0: data do not have skewness")
    print("pvalsk")
    print(pvalsk)
    pskneg <- c(numeric(nvars))
    for (p in 1:nvars) pskneg[p] <- pnorm(z1[p])
    print("H0: data do not have negative skewness")
    print("pskneg")
    print(pskneg)
    pskpos <- 1 - pskneg
    print("H0: data do not have positive skewness")
    print("pskpos")
    print(pskpos)
    pvalk <- c(numeric(nvars))
    for (p in 1:nvars) pvalk[p] <- pnorm(z2[p], lower.tail = FALSE)
    for (p in 1:nvars) pvalk[p] <- 2 * pvalk[p]
    for (p in 1:nvars) if (pvalk[p] > 1) 
        pvalk[p] <- 2 - pvalk[p]
    print("H0: data do not have kurtosis")
    print("pvalk")
    print(pvalk)
    pkneg <- c(numeric(nvars))
    for (p in 1:nvars) pkneg[p] <- pnorm(z2[p])
    print("H0: data do not have negative kurtosis")
    print("pkneg")
    print(pkneg)
    pkpos <- 1 - pkneg
    print("H0: data do not have positive kurtosis")
    print("pkpos")
    print(pkpos)
    z1 <- matrix(z1, nrow = 1)
    z2 <- matrix(z2, nrow = 1)
    Ep <- z1 %*% t(z1) + z2 %*% t(z2)
    dof <- 2 * nvars
    sig.Ep <- 1 - pchisq(Ep, dof)
    print("H0: data are normally distributed")
    print("Ep")
    print(Ep)
    print("dof")
    print(dof)
    print("sig.Ep")
    print(sig.Ep)
}
