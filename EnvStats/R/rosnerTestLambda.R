rosnerTestLambda <-
function (n, k = 10, alpha = 0.05) 
{
    if (!is.numeric(n) || !is.finite(n) || !all(n == round(n)) || 
        any(n < 3)) 
        stop("All values of 'n' must be positive integers greater than 2")
    if (!is.numeric(k) || !is.finite(k) || !all(k == round(k)) || 
        any(k < 1) || any(k > (n - 2))) 
        stop("All values of 'k' must be positive integers less than or equal to n-2")
    if (!is.numeric(alpha) || !is.finite(alpha) || any(alpha <= 
        0) || any(alpha >= 1)) 
        stop("All values of 'alpha' must be greater than 0 and less than 1")
    arg.mat <- cbind.no.warn(n = as.vector(n), k = as.vector(k), 
        alpha = as.vector(alpha))
    for (i in c("n", "k", "alpha")) assign(i, arg.mat[, i])
    l <- k - 1
    p <- 1 - ((alpha/2)/(n - l))
    t.crit <- qt(p = p, df = n - l - 2)
    lambda <- t.crit * (n - l - 1)/sqrt((n - l - 2 + t.crit^2) * 
        (n - l))
    names(lambda) <- NULL
    lambda
}
