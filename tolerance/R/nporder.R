np.order <- function (m, alpha = 0.05, P = 0.99, indices = FALSE) 
{
    f <- function(n, m, P, alpha) pbeta(1 - P, m, n - m + 1) - 
        (1 - alpha)
    n <- try(suppressWarnings(ceiling(uniroot(f, interval = c(m, 
        1e+05 * m), m = m, alpha = alpha, P = P)$root)), silent = TRUE)
    if (class(n) == "try-error") 
        stop("Interval to search for the root is too wide.", 
            call. = FALSE)
    if (indices == FALSE | (m < 2)) {
        print(n)
    }
    else {
        temp <- list()
        temp[[1]] <- n
        ind <- NULL
        for (r in 1:(m - 1)) {
            ind <- rbind(ind, c(r, n - (m - r) + 1))
        }
        temp[[2]] <- ind
        names(temp) <- c("Sample Size", "Order Statistics")
        temp
    }
}
