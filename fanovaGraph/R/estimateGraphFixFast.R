estimateGraphFixFast <- function(f.mat, d, q, q.arg, n.mc, n.fast, print.loop.index, ...) {
    # FAST frequencies
    w <- c(11, 35)
    if (n.fast < 2 * 6 * max(w)) 
        stop("n.fast too small to guarantee positive values")
    fast <- function(f, ...) {
        Y <- drop(f(X, ...))
        Var <- 2 * sum(colMeans(Y * sinsop)^2 + colMeans(Y * cossop)^2)
        D1 <- 2 * sum(colMeans(Y * sinsop[, harm1])^2 + colMeans(Y * 
            cossop[, harm1])^2)
        D2 <- 2 * sum(colMeans(Y * sinsop[, harm2])^2 + colMeans(Y * 
            cossop[, harm2])^2)
        return(Var - D1 - D2)
    }
    ## the fixed values, only d-2 points needed, but programming much
    #   shorter with d
    SampleFixed <- matrix(runif(n.mc * d), ncol = d)
    for (j in 1:d) SampleFixed[, j] <- do.call(q[j], c(list(p = SampleFixed[, 
        j]), q.arg[[j]]))
    # function f.mat but depending only on Xjk, other variables fixed
    #   by xfixed
    fjk.mat <- function(Xjk, jk, xfixed, ...) {
        X <- rep(1, nrow(Xjk)) %*% t(xfixed)
        X[, jk[1]] <- Xjk[, 1]
        X[, jk[2]] <- Xjk[, 2]
        return(f.mat(X, ...))
    }
    JK <- t(combn(1:d, 2))
    DintJK <- numeric(dim(JK)[1])
    # Design matrix for fast:
    s <- seq(-pi, pi, len = n.fast)
    p <- 1:210
    sop <- s %o% p
    sinsop <- sin(sop)
    cossop <- cos(sop)
    harm1 <- 1:6 * w[1]
    harm2 <- 1:6 * w[2]
    X01 <- matrix(, nrow = n.fast, ncol = 2)
    for (i in 1:2)
      X01[, i] <- 1/2 + 1/pi * asin(sin(w[i] * s))
    X <- matrix(, nrow = n.fast, ncol = 2)
    for (j in 1:nrow(JK)) {
        for (i in 1:2) {
           X[, i] <- do.call(q[JK[j,i]], c(list(p = X01[, i]), q.arg[[JK[j,i]]]))
        }
      Dint <- numeric(n.mc)
        if(print.loop.index) cat("index = ", JK[j, ],"\n") # collapse=""
        for (m in (1:n.mc)) {
            Dint[m] <- fast(fjk.mat, jk = JK[j, ], xfixed = SampleFixed[m, ], ...)
        }
        DintJK[j] <- mean(Dint)
    }
    inter <- paste("X",JK[,1],"*","X",JK[,2], sep="")
    totalInt <- as.matrix(DintJK)
    rownames(totalInt) <- inter
    return(totalInt)
} 
