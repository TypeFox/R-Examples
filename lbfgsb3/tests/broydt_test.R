rm(list = ls())
library(lbfgsb3)

broydt.f <- function(x) {
    n <- length(x)
    res <- rep(NA, n)
    res[1] <- ((3 - 0.5 * x[1]) * x[1]) - 2 * x[2] + 1
    tnm1 <- 2:(n - 1)
    res[tnm1] <- ((3 - 0.5 * x[tnm1]) * x[tnm1]) - x[tnm1 - 1] - 
        2 * x[tnm1 + 1] + 1
    res[n] <- ((3 - 0.5 * x[n]) * x[n]) - x[n - 1] + 1
    sum(res * res)
}

broydt.g <- function(x) {
    n <- length(x)
    gg <- rep(NA, n)  # gradient set to NA to start with
    gg[1] <- -2 + 2 * x[1] + 4 * x[3] + (6 - 2 * x[1]) * (1 - 
        2 * x[2] + x[1] * (3 - 0.5 * x[1])) - 2 * x[2] * (3 - 
        0.5 * x[2])
    gg[2] <- -6 + 4 * x[4] + 10 * x[2] + (6 - 2 * x[2]) * (1 - 
        x[1] - 2 * x[3] + x[2] * (3 - 0.5 * x[2])) - 4 * x[1] * 
        (3 - 0.5 * x[1]) - 2 * x[3] * (3 - 0.5 * x[3])
    tnm2 <- 3:(n - 2)
    gg[tnm2] <- -6 + 4 * x[tnm2 - 2] + 4 * x[tnm2 + 2] + 10 * 
        x[tnm2] + (6 - 2 * x[tnm2]) * (1 - x[tnm2 - 1] - 2 * 
        x[tnm2 + 1] + x[tnm2] * (3 - 0.5 * x[tnm2])) - 4 * x[tnm2 - 
        1] * (3 - 0.5 * x[tnm2 - 1]) - 2 * x[tnm2 + 1] * (3 - 
        0.5 * x[tnm2 + 1])
    gg[n - 1] <- -6 + 4 * x[n - 3] + 10 * x[n - 1] + (6 - 2 * 
        x[n - 1]) * (1 - x[n - 2] - 2 * x[n] + x[n - 1] * (3 - 
        0.5 * x[n - 1])) - 4 * x[n - 2] * (3 - 0.5 * x[n - 2]) - 
        2 * x[n] * (3 - 0.5 * x[n])
    
    gg[n] <- -4 + 4 * x[n - 2] + 8 * x[n] + (6 - 2 * x[n]) * 
        (1 - x[n - 1] + x[n] * (3 - 0.5 * x[n])) - 4 * x[n - 
        1] * (3 - 0.5 * x[n - 1])
    return(gg)
}

ni <- c(10, 100, 400)

for (ii in ni) {
    n <- ii
    cat("n=", n, "\n")
    x0 <- rep(pi, n)
    ut <- system.time(ans <- lbfgsb3(x0, broydt.f, broydt.g, control = list(trace = 1)))[1]
    cat("unconstrained n=",n," with gradient\n")
    cat("time =",ut,"\n")
    print(ans)

    ct <- system.time(ans <- lbfgsb3(x0, broydt.f))[1]
    cat("unconstrained n=",n," no analytic gradient\n")
    cat("time =",ct,"\n")
    print(ans)
    
    
    lower <- rep(1, n)
    upper <- rep(Inf, n)

    x0 <- rep(pi, n)
    t1 <- system.time(ans <- lbfgsb3(x0, broydt.f, broydt.g, lower = lower, 
        upper = upper))[1]
    cat("constrained n= ",n," no analytic gradient\n")
    cat("time =",t1,"\n")
    print(ans)

}

