weightLine <- function(p, w, col0, lwd0, lty0 = 1, type = c("pval", "empirical")[1]){
    n <- length(p)
    k <- length(w)
    ps <- pPool(p)
    ws <- w
    lower <- c(ps[2:k], 0)  
    odd1 <- n %% 2
    even1 <- (n + 1) %% 2
    upper <- c(ps[1:(k - 1)], odd1 * rev(sort(p))[n - 1] + even1 * rev(sort(p))[n])

    if (type == "empirical"){
        lower <- (k - 1):0 / k
        upper <- (k:1) / k
    }

    segments(lower, ws, upper, ws, col = col0, lwd = lwd0, lty = lty0)
    segments(upper[k:2], ws[(k - 1):1], upper[k:2], ws[k:2], col = col0, lwd = lwd0, lty = lty0)
}
