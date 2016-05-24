snut <- function(Y, X,
                 cluster = rep(1, length(Y)),
                 offset = rep(0, length(Y))
                 ){
    ## X may not have a constant column!
    n <- length(Y)
    if (is.matrix(X)){
        if (NROW(X) != n) stop("Wrong dimension of X")
        q <- NCOL(X)
    }else{
        if (length(X) != n) stop("Wrong length of X")
        q <- 1
        X <- matrix(X, ncol = 1)
    }

    cluster <- as.vector(unclass(factor(cluster)))

    in.here <- function(y) (sum(y) > 0.5) && (sum(y) < (length(y) - 0.5))
    
    here <- cluster %in% which(tapply(Y, cluster, in.here))
    Y <- Y[here]
    X <- X[here, , drop = FALSE]
    cluster <- cluster[here]
    offset <- offset[here]
    
    p <- length(unique(cluster))
    pq <- p + q
    beta <- numeric(pq)

    cluster <- as.vector(unclass(factor(cluster)))
    
    fun <- function(beta){
        ## Minus The log likelihood:
        lin <- offset + beta[cluster] + X %*% beta[(p+1):pq]
        -sum(Y * lin - log(1 + exp(lin)))
    }

    gra <- function(beta){
        # Minus the gradient
        ret <- numeric(length(beta))
        lin <- offset + beta[cluster] + X %*% beta[(p+1):pq]
        P <- exp(lin)
        P <- P / (1 + P)
        ymP <- Y - P # is an nx1 matrix
        ret[1:p] <- tapply(ymP, cluster, sum)
        ret[(p+1):pq] <- t(ymP) %*% X
        ##for (s in 1:q) ret[s + p] <- sum(X[, s] * ymP)
        -ret
    }

    res <- optim(beta, fun, gra, method = "BFGS", control = list(trace = TRUE))
    res
}
