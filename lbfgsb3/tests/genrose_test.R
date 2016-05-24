## Optimization test function GENROSE
## ?? refs (put in .doc??)
rm(list = ls())
library(lbfgsb3)

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



cat("\n\n Unconstrained test\n")
nn <- 100
xx <- rep(3, nn)
lo <- -Inf
up <- Inf
t100u <- system.time(ans100u <- lbfgsb3(xx, genrose.f, 
    genrose.g, gs = 10))[1]
cat("final fn value =", ans100u$f," after ",ans100u$info$isave[34], "\n")
cat("time = ", t100u, "\n")
t100uo <- system.time(ao100u <- optim(xx, genrose.f, 
    genrose.g, method = "L-BFGS-B", gs = 10, control=list(trace=1)))[1]
cat("final fn value =", ao100u$value," after",ao100u$counts[1],ao100u$counts[2], "\n")
cat("time = ", t100uo, "\n")

t100un <- system.time(ans100un <- lbfgsb3(xx, genrose.f, 
    gr = NULL, gs = 10))[1]
cat("final fn value =", ans100un$f, " after ",ans100un$info$isave[34],"\n")
cat("time = ", t100un, "\n")
t100uon <- system.time(ao100un <- optim(xx, genrose.f, 
    gr = NULL, method = "L-BFGS-B", gs = 10, control=list(trace=1)))[1]
cat("final fn value =", ao100un$value," after",ao100un$counts[1],ao100un$counts[2],  "\n")
cat("time = ", t100uon, "\n")
