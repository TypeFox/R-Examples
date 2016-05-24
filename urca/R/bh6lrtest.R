##
## bh6lrtest
##
bh6lrtest <- function (z, H, r, r1, conv.val=0.0001, max.iter=50) 
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
    s <- ncol(H)
    r2 <- r - r1
    lambda <- z@lambda
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
    Dtemp <- chol(t(H)%*%SKK%*%H, pivot = TRUE)
    pivot <- attr(Dtemp, "pivot")
    oo <- order(pivot)
    D <- t(Dtemp[, oo])
    Dinv <- solve(D)
    valeigen <- eigen(Dinv %*% t(H) %*% SK0 %*% solve(S00) %*% S0K %*% H %*% t(Dinv))
    beta1 <- H%*%valeigen$vectors[,1:r1]
    i <- 0
    last <- 1
    diff <- 1
    while(diff > conv.val){
      S00.b1 <- S00 - S0K%*%beta1%*%solve(t(beta1)%*%SKK%*%beta1)%*%t(beta1)%*%SK0
      S0K.b1 <- S0K - S0K%*%beta1%*%solve(t(beta1)%*%SKK%*%beta1)%*%t(beta1)%*%SKK
      SK0.b1 <- SK0 - SKK%*%beta1%*%solve(t(beta1)%*%SKK%*%beta1)%*%t(beta1)%*%SK0
      SKK.b1 <- SKK - SKK%*%beta1%*%solve(t(beta1)%*%SKK%*%beta1)%*%t(beta1)%*%SKK
      valeigen <- eigen(SKK.b1)
      C <- valeigen$vectors[ ,1:(P-r1)]%*%diag(1/sqrt(valeigen$values[1:(P-r1)]))
      valeigen <- eigen(t(C)%*%SK0.b1%*%solve(S00.b1)%*%S0K.b1%*%C)
      lambda.res <- valeigen$values
      diff <- t(lambda.res-last)%*%(lambda.res-last)
      last <- lambda.res
      beta2 <- C%*%valeigen$vectors[,1:r2]
      S00.b2 <- S00 - S0K%*%beta2%*%solve(t(beta2)%*%SKK%*%beta2)%*%t(beta2)%*%SK0
      S0K.b2 <- S0K - S0K%*%beta2%*%solve(t(beta2)%*%SKK%*%beta2)%*%t(beta2)%*%SKK
      SK0.b2 <- SK0 - SKK%*%beta2%*%solve(t(beta2)%*%SKK%*%beta2)%*%t(beta2)%*%SK0
      SKK.b2 <- SKK - SKK%*%beta2%*%solve(t(beta2)%*%SKK%*%beta2)%*%t(beta2)%*%SKK
      valeigen <- eigen(t(H)%*%SKK.b2%*%H)
      C <- valeigen$vectors[ ,1:s]%*%diag(1/sqrt(valeigen$values[1:s]))
      valeigen <- eigen(t(C)%*%t(H)%*%SK0.b2%*%solve(S00.b2)%*%S0K.b2%*%H%*%C)
      beta1 <- H%*%valeigen$vectors[,1:r1]
      i <- i + 1
      if(i>max.iter){
        warning("\nNo convergence, used last iterations values.\n")
        break
      }
    }
    Vorg <- cbind(beta1, beta2)
    V <- Vorg
    idx <- ncol(V)
    V <- sapply(1:idx, function(x) V[, x]/V[1, x])
    W <- S0K %*% V %*% solve(t(V) %*% SKK %*% V)
    PI <- W %*% t(V)
    DELTA <- S00 - S0K %*% V %*% solve(t(V) %*% SKK %*% V) %*% t(V) %*% SK0
    GAMMA <- M01 %*% M11inv - PI %*% MK1 %*% M11inv
    Dtemp <- chol(t(beta1)%*%SKK%*%beta1, pivot = TRUE)
    pivot <- attr(Dtemp, "pivot")
    oo <- order(pivot)
    D <- t(Dtemp[, oo])
    Dinv <- solve(D)
    valeigen <- eigen(Dinv %*% t(beta1) %*% SK0 %*% solve(S00) %*% S0K %*% beta1 %*% t(Dinv))
    rho <- valeigen$values
    teststat <- N*(sum(log(1-rho[1:r1])) + sum(log(1-lambda.res[1:r2])) - sum(log(1-lambda[1:r])))
    df <- (P - s - r2)*r1
    pval <- c(1 - pchisq(teststat, df), df)
    new("cajo.test", Z0 = z@Z0, Z1 = z@Z1, ZK = z@ZK, ecdet = z@ecdet, H = H, A = NULL, B = NULL, type = type, teststat = teststat, pval = pval, lambda = lambda.res, Vorg = Vorg, V = V, W = W, PI = PI, DELTA = DELTA, DELTA.bb = NULL, DELTA.ab = NULL, DELTA.aa.b = NULL, GAMMA = GAMMA, test.name = "Johansen-Procedure")
}
