library(testthat)

context("wkm")

## TODO: test extreme cases (no censoring, all obs. censored, ...)

## TODO: test with data from Kalbfleisch/Prentice paper

test_that("single stratum wkm estimate reduces to Kaplan-Meier estimator", {
    test.data <- function(n) {
        Z <- rep(0, n)
        T <- rexp(n, 0.1)
        C <- runif(n, 0, 10)
        X <- pmin(T, C)
        D <- T <= C
        V <- runif(n, 0, X/4) 
        
        data.frame(Y=X, D=D, W=Z, V=V)
  }

    data <- test.data(100)    
    times <- seq(0, 5, length.out=10) ##c(1.5, 3)

    data <- data[order(data$Y),]
  
    S <- exp(-times/10)

    fs <- survfit(Surv(V, Y, D) ~ 1, data=data)
    fit <- wkm(times, formula=Surv(V, Y, D) ~ strata(W), data=data)
    
    f <- approxfun(fs$time, fs$surv, method="constant", yleft=1, rule=2, f=0)
    g <- approxfun(fs$time, 100 * (fs$surv * fs$std.err)^2, method="constant", yleft=0, rule=2, f=0)

    expect_true(all.equal(f(times), fit$S))
    expect_true(all.equal(g(times), fit$V))
})

test_that("covariance matrix is symmetric", {
    X <- sample(c(0,1), 100, replace=TRUE)
    T <- numeric(100)
    T[X == 0] <- rexp(sum(X == 0), 0.1)
    T[X == 1] <- rexp(sum(X == 1), 0.5)
    
    C <- runif(100, 0, 10)
    
    Y <- pmin(T, C)
    D <- T <= C
    V <- runif(100, 0, Y/2)
  
    data <- data.frame(Y=Y, D=D, W=X, V=V)
    
    expect_true(isSymmetric(wkm(sort(Y), formula=Surv(V, Y, D) ~ strata(W), data=data, param=list(cov=TRUE, var=TRUE, alpha=1, left.limit=FALSE))$COV))
})

test_that("variance and covariance calculations match", {
    T <- rexp(200)
    C <- rexp(200)
    D <- T <= C
    Y <- pmin(T, C)

    Z <- rbinom(200, 1, 0.5)
    X <- rbinom(200, 1, 0.5)

    data <- data.frame(Y=Y, D=D, Trt=Z, W=factor(X))

    times <- sort(data$Y[data$D == 1])
    
    fit <- wkm(times, formula=Surv(Y, D) ~ strata(W), data=data, param=list(cov=TRUE, var=FALSE, alpha=1, left.limit=FALSE))
    fit2 <- wkm(times, formula=Surv(Y, D) ~ strata(W), data=data, param=list(cov=FALSE, var=TRUE, alpha=1, left.limit=FALSE))
    
    expect_true(all.equal(diag(fit$COV), fit2$V))
})
