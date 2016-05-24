##
## bh5lrtest
##
bh5lrtest <- function (z, H, r) 
{
    if (!(class(z) == "ca.jo")) {
        stop("\nPlease, provide object of class 'ca.jo' as 'z'.\n")
    }
    if (r >= z@P || r < 1) {
        stop("\nCount of cointegrating relationships is out of allowable range.\n")
    }
    if (z@ecdet == "none") {
        P <- z@P
    } else {
        P <- z@P + 1
    }
    r <- as.integer(r)
    H <- as.matrix(H)
    if (!(nrow(H) == P)) {
        stop("\nRow number of 'H' is unequal to VAR order.\n")
    }
    if(ncol(H) > r - 1){
      stop("\nToo many columns in H for provided r.\n")
    }
    r1 <- ncol(H)
    r2 <- r - r1
    type <- "Estimation and testing under partly known beta"
    N <- nrow(z@Z0)
    M00 <- crossprod(z@Z0)/N
    M11 <- crossprod(z@Z1)/N
    MKK <- crossprod(z@ZK)/N
    M01 <- crossprod(z@Z0, z@Z1)/N
    M0K <- crossprod(z@Z0, z@ZK)/N
    MK0 <- crossprod(z@ZK, z@Z0)/N
    M10 <- crossprod(z@Z1, z@Z0)/N
    M1K <- crossprod(z@Z1, z@ZK)/N
    MK1 <- crossprod(z@ZK, z@Z1)/N
    M11inv <- solve(M11)
    S00 <- M00 - M01 %*% M11inv %*% M10
    S0K <- M0K - M01 %*% M11inv %*% M1K
    SK0 <- MK0 - MK1 %*% M11inv %*% M10
    SKK <- MKK - MK1 %*% M11inv %*% M1K
    S00.h <- S00 - S0K%*%H%*%solve(t(H)%*%SKK%*%H)%*%t(H)%*%SK0
    S0K.h <- S0K - S0K%*%H%*%solve(t(H)%*%SKK%*%H)%*%t(H)%*%SKK
    SK0.h <- SK0 - SKK%*%H%*%solve(t(H)%*%SKK%*%H)%*%t(H)%*%SK0
    SKK.h <- SKK - SKK%*%H%*%solve(t(H)%*%SKK%*%H)%*%t(H)%*%SKK
    valeigen <- eigen(SKK.h)
    C <- valeigen$vectors[ ,1:(P-r1)]%*%diag(1/sqrt(valeigen$values[1:(P-r1)]))
    valeigen <- eigen(t(C)%*%SK0.h%*%solve(S00.h)%*%S0K.h%*%C)
    PSI <- C%*%valeigen$vectors[,1:r2]
    Dtemp <- chol(t(H)%*%SKK%*%H, pivot = TRUE)
    pivot <- attr(Dtemp, "pivot")
    oo <- order(pivot)
    D <- t(Dtemp[, oo])
    Dinv <- solve(D)
    rho <- eigen(Dinv %*% t(H) %*% SK0 %*% solve(S00) %*% S0K %*% H %*% t(Dinv))
    Vorg <- cbind(H, PSI)
    idx <- ncol(PSI)
    PSI <- sapply(1:idx, function(x) PSI[, x]/PSI[1, x])
    V <- cbind(H, PSI)
    W <- S0K %*% V %*% solve(t(V) %*% SKK %*% V)
    PI <- W %*% t(V)
    DELTA <- S00 - S0K %*% V %*% solve(t(V) %*% SKK %*% V) %*% t(V) %*% SK0
    GAMMA <- M01 %*% M11inv - PI %*% MK1 %*% M11inv
    lambda.res <- valeigen$values
    lambda <- z@lambda
    teststat <- N *(sum(log(1 - rho$values[1:r1])) + sum(log(1-lambda.res[1:r2])) - sum(log(1 - lambda[1:r])))
    df <- (P - r)*r1
    pval <- c(1 - pchisq(teststat, df), df)
    new("cajo.test", Z0 = z@Z0, Z1 = z@Z1, ZK = z@ZK, ecdet = z@ecdet, H = H, A = NULL, B = NULL, type = type, teststat = teststat, pval = pval, lambda = lambda.res, Vorg = Vorg, V = V, W = W, PI = PI, DELTA = DELTA, DELTA.bb = NULL, DELTA.ab = NULL, DELTA.aa.b = NULL, GAMMA = GAMMA, test.name = "Johansen-Procedure")
}
