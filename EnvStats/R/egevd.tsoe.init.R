egevd.tsoe.init <-
function (x, p) 
{
    C <- -log(p)
    A13 <- C[1]/C[3]
    A23 <- C[2]/C[3]
    D <- (x[2] - x[3])/(x[1] - x[3])
    if (D < (log(A23)/log(A13))) {
        A21 <- C[2]/C[1]
        lower <- sqrt(.Machine$double.eps)
        upper <- log(D)/log(A21)
    }
    else {
        lower <- log(1 - D)/log(A23)
        upper <- -sqrt(.Machine$double.eps)
    }
    k <- uniroot(egevd.tsoe.init.h, lower = lower, upper = upper, 
        A13 = A13, A23 = A23, D = D)$root
    alpha <- (k * (x[1] - x[3]))/(C[3]^k - C[1]^k)
    lambda <- x[1] - (alpha * (1 - C[1]^k))/k
    c(lambda, alpha, k)
}
