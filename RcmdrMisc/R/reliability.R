# last modified 2014-08-04 by J. Fox

reliability <- function(S){
    reliab <- function(S, R){
        k <- dim(S)[1]
        ones <- rep(1, k)
        v <- as.vector(ones %*% S %*% ones)
        alpha <- (k/(k - 1)) * (1 - (1/v)*sum(diag(S)))
        rbar <- mean(R[lower.tri(R)])
        std.alpha <- k*rbar/(1 + (k - 1)*rbar)
        c(alpha=alpha, std.alpha=std.alpha)
    }
    result <- list()
    if ((!is.numeric(S)) || !is.matrix(S) || (nrow(S) != ncol(S))
        || any(abs(S - t(S)) > max(abs(S))*1e-10) || nrow(S) < 2)
        stop("argument must be a square, symmetric, numeric covariance matrix")
    k <- dim(S)[1]
    s <- sqrt(diag(S))
    R <- S/(s %o% s)
    rel <- reliab(S, R)
    result$alpha <- rel[1]
    result$st.alpha <- rel[2]
    if (k < 3) {
        warning("there are fewer than 3 items in the scale")
        return(invisible(NULL))
    }
    rel <- matrix(0, k, 3)
    for (i in 1:k) {
        rel[i, c(1,2)] <- reliab(S[-i, -i], R[-i, -i])
        a <- rep(0, k)
        b <- rep(1, k)
        a[i] <- 1
        b[i] <- 0
        cov <- a %*% S %*% b
        var <- b %*% S %*% b
        rel[i, 3] <- cov/(sqrt(var * S[i,i]))
    }
    rownames(rel) <- rownames(S)
    colnames(rel) <- c("Alpha", "Std.Alpha", "r(item, total)")
    result$rel.matrix <- rel
    class(result) <- "reliability"
    result
}

print.reliability <- function(x, digits=4, ...){
    cat(paste("Alpha reliability = ", round(x$alpha, digits), "\n"))
    cat(paste("Standardized alpha = ", round(x$st.alpha, digits), "\n"))
    cat("\nReliability deleting each item in turn:\n")
    print(round(x$rel.matrix, digits))
    invisible(x)
}
