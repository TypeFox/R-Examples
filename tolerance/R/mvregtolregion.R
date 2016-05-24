mvregtol.region <- function (y, x, new.x = NULL, int = TRUE, alpha = 0.05, 
    P = 0.99, B = 1000) 
{
    temp.x <- colnames(x)
    temp.y <- colnames(y)
    x <- as.matrix(x)
    q <- ncol(y)
    n <- nrow(y)
    m <- ncol(x)
    new.x <- matrix(new.x, ncol = m)
    f.m <- n - m - 1
    P.n <- matrix(1, n, n)
    P.n <- diag(1, n) - (P.n/n)
    x.all <- rbind(x, new.x)
    N <- nrow(x.all)
    x.bar <- apply(x, 2, mean)
    A <- try(solve(t(x) %*% P.n %*% x), silent = TRUE)
    if (class(A) == "try-error") 
        stop(paste("The design matrix is not of full rank!", 
            "\n"))
    if (int) {
        X <- cbind(1, x)
        X.all <- cbind(1, x.all)
    }
    else {
        X <- x
        X.all <- x.all
    }
    B.hat <- solve(t(X) %*% X) %*% t(X) %*% y
    y.hat <- X.all %*% B.hat
    d.2 <- (1/n) + sapply(1:N, function(i) rbind(x.all[i, ] - 
        x.bar) %*% A %*% cbind(x.all[i, ] - x.bar))
    H.2 <- lapply(1:B, function(i) matrix(rchisq(N * q, df = 1), 
        ncol = q) * d.2)
    L <- t(sapply(1:B, function(i) eigen(rwishart(f.m, q))$values))
    c1 <- sapply(1:B, function(i) apply((1 + H.2[[i]]^2)/L[[i]], 
        1, sum))
    c2 <- sapply(1:B, function(i) apply((1 + 2 * (H.2[[i]])^2)/(L[[i]])^2, 
        1, sum))
    c3 <- sapply(1:B, function(i) apply((1 + 3 * (H.2[[i]])^2)/(L[[i]])^3, 
        1, sum))
    a <- (c2^3)/(c3^2)
    T.all <- lapply(1:N, function(i) f.m * ((sqrt(c2[i, ])/a[i, 
        ]) * (qchisq(P, a[i, ]) - a[i, ]) + c1[i, ]))
    k <- sapply(T.all, quantile, 1 - alpha)
    cat("These are ", (1 - alpha) * 100, "%/", P * 100, "% tolerance factors.", 
        sep = "", fill = TRUE)
    tol <- cbind(k, data.frame(y.hat), data.frame(x.all))
    if (is.null(temp.x)) 
        temp.x <- paste("x", 1:m, sep = "")
    if (is.null(temp.y)) {
        temp.y <- paste("y", 1:q, ".hat", sep = "")
    }
    else temp.y <- paste(temp.y, ".hat", sep = "")
    tol.col <- c("k.factor", temp.y, temp.x)
    colnames(tol) <- tol.col
    tol
}
