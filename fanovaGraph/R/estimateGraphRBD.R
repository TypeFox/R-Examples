estimateGraphRBD <- function(f.mat, d, q, q.arg, L, M, print.loop.index, ...) {
    if (L > 4000) 
        stop("L < 4000 required")
    N <- 2 * (M * d + L)
    w <- d:1
    s <- -pi + 2 * pi/N * (1:N)
    X <- matrix(, N, d)  # Design X
    for (i in 1:d) {
        X[, i] <- 1/2 + 1/pi * asin(sin(w[i] * s))
        X[, i] <- do.call(q[i], c(list(p = X[, i]), q.arg[[i]]))
    }
    ps <- (d * M + 1):(N/2)  # values of p
    sop <- s %o% ps  # matrix s * p
    sin_sop <- sin(sop)
    cos_sop <- cos(sop)
    
    RBDt <- function(f.mat, d, order, L, M, ...) {
        JK <- combn(1:d, order)  # factor combinations
        DTi <- numeric(ncol(JK))
        for (i in 1:ncol(JK)) {
            # for all factor combinations
            if(print.loop.index) cat("index = ", JK[,i],"\n") 
            Xs <- X
            o <- sample(1:N)  # X sampled at i
            Xs[, JK[, i]] <- Xs[o, JK[, i]]
            Y <- drop(f.mat(Xs, ...))
            Amp <- colMeans(Y * cos_sop)^2 + colMeans(Y * sin_sop)^2
            DTi[i] <- N/L * sum(Amp)
        }
        rbind(JK, DTi)
    }
    STij <- RBDt(f.mat, d, order = 2, L = L, M = M, ...)
    STi <- RBDt(f.mat, d, order = 1, L = L, M = M, ...)
    
    totalInt <- STij  # just initializing
    for (i in 1:ncol(totalInt)) totalInt[3, i] <- sum(STi[2, totalInt[1:2, 
        i]]) - STij[3, i]
    inter <- paste("X",totalInt[1,],"*","X",totalInt[2,], sep="")
    totalInt <- as.matrix(totalInt[3,])
    rownames(totalInt) <- inter
    return(totalInt)
} 