acc.samp <- function (n, N, alpha = 0.05, P = 0.99, AQL = 0.01, RQL = 0.02) 
{
    if (RQL - AQL < 4e-08) {
        stop(paste("RQL must be greater than AQL!"))
    }
    D <- N - (N * P)
    m.h <- D
    n.h <- N - D
    ff <- function(k, c, m, n) (alpha) - phyper(c, m, n, k)
    c <- try(floor(uniroot(ff, interval = c(0, ceiling(D)), k = n, m = m.h, 
        n = n.h)$root), silent = TRUE)
    if (class(c) == "try-error") 
        c <- 0
    if (phyper(c, m = m.h, n = n.h, k = n) > alpha) 
        c <- max(c - 1, 0)
    prob.rej.good <- 1 - phyper(c, m=floor(AQL * N), n=N - floor(AQL * 
        N), k=n)
    prob.rej.bad <- 1 - phyper(c, m=floor(RQL * N), n=(N - floor(RQL * 
        N)), k = n)
    temp <- c(round(c, 0), round(N, 0), 1 - round(alpha, 
        4), round(P,4), round(AQL, 4), round(RQL, 4), round(n, 0), round(prob.rej.good, 
        4), round((1 - prob.rej.bad), 4))
    temp <- as.matrix(temp)
    rownames(temp) <- c("acceptance.limit", "lot.size", "confidence", "P", "AQL", 
         "RQL", "sample.size", "prod.risk", "cons.risk")
    colnames(temp) <- ""
    if (temp[8, ] > alpha) 
        cat("Warning: Desired confidence level not attained!", 
            "\n")
    if (n == N)
        cat("Sample size and lot size are the same.  This is just a census!",
            "\n")
    temp
}
