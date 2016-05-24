## Optimization test function GENROSE
## ?? refs (put in .doc??)
rm(list = ls())
library(Rcgmin)

genrose.f <- function(x, gs = NULL) {
    # objective function
    ## One generalization of the Rosenbrock banana valley
    #   function (n parameters)
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    fval <- 1 + sum(gs * (x[1:(n - 1)]^2 - x[2:n])^2 + (x[2:n] - 
        1)^2)
    return(fval)
}

genrose.g0 <- function(x, gs = NULL) {
    # gradient for genrosef, genrose.f / genrosep.f
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    gg <- as.vector(rep(0, n))
    for (i in 2:n) {
        z1 <- x[i] - x[i - 1] * x[i - 1]
        z2 <- 1 - x[i]
        gg[i] <- 2 * (gs * z1 - z2)
        gg[i - 1] <- gg[i - 1] - 4 * gs * x[i - 1] * z1
    }
    return(gg)
}

genrose.h <- function(x, gs = NULL) {
    ## compute Hessian
    if (is.null(gs)) {
        gs = 100
    }
    n <- length(x)
    hh <- matrix(rep(0, n * n), n, n)
    for (i in 2:n) {
        z1 <- x[i] - x[i - 1] * x[i - 1]
        z2 <- 1 - x[i]
        hh[i, i] <- hh[i, i] + 2 * (gs + 1)
        hh[i - 1, i - 1] <- hh[i - 1, i - 1] - 4 * gs * z1 - 
            4 * gs * x[i - 1] * (-2 * x[i - 1])
        hh[i, i - 1] <- hh[i, i - 1] - 4 * gs * x[i - 1]
        hh[i - 1, i] <- hh[i - 1, i] - 4 * gs * x[i - 1]
    }
    return(hh)
}

genrose.doc <- function() {
    ## documentation for genrose
    cat("One generalization of the Rosenbrock banana valley function (n parameters)\n")
    ## How should we do the documentation output?
}

genrose.g <- function(x, gs = NULL) {
    # vectorized gradient for genrose.f
    # Ravi Varadhan 2009-04-03
    n <- length(x)
    if (is.null(gs)) {
        gs = 100
    }
    gg <- as.vector(rep(0, n))
    tn <- 2:n
    tn1 <- tn - 1
    z1 <- x[tn] - x[tn1]^2
    z2 <- 1 - x[tn]
    gg[tn] <- 2 * (gs * z1 - z2)
    gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1
    gg
}

##  Some timing results
xx <- rep(2, 10000)
##k100v<-system.time(gra<-genrose.g(xx,gs=100.0))
##k100v0<-system.time(gra<-genrose.g0(xx,gs=100.0))
##  k100v
##    user  system elapsed
##    0.068   0.016   0.084
##  k100v0
##     user  system elapsed
##    3.188   0.004   3.198
##


t10k <- system.time(ans <- Rcgmin(xx, genrose.f, genrose.g0, 
    control = list(trace = 1), gs = 100))[1]
cat("final fn value =", ans$value, "\n")
cat("time = ", t10k, "\n")


cat("\n\n Masked test \n")
xx <- rep(2, 1000)
bdmsk <- c(rep(0, 10), rep(1, 990))
ll <- rep(-20, 1000)
uu <- -ll
t1000m <- system.time(ans <- Rcgmin(xx, genrose.f, 
    genrose.g0, lower = ll, upper = uu, bdmsk = bdmsk, control = list(trace = 1), 
    gs = 100))[1]
cat("final fn value =", ans$value, "\n")
cat("time = ", t1000m, "\n")
