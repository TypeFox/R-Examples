source("helper-diversitree.R")

context("Matrix exponentials")

set.seed(1)
n <- 64L
Q <- diversitree:::mkn.Q(runif(n*(n-1)))
v <- runif(n)
v <- v / sum(v)
t <- 0.02

## 1. Plain matrix exponential
E.f <- diversitree:::expm.dense(Q, t)
E.e <- expm(Q*t)
expect_that(E.f, equals(E.e))

## 2. exp(Qt)v
ans.f <- diversitree:::expmv.expokit.dense(Q, t, v)
ans.e <- c(expm(Q*t) %*% v)
expect_that(ans.f, equals(ans.e))

## 3. sparse exp(Qt)v
set.seed(1)
i <- sort(sample(n*(n-1), 2*n))
p <- rep(0.0, n*(n-1))
p[i] <- runif(2*n)
Q2 <- diversitree:::mkn.Q(p)

ans.e <- c(expm(Q2*t) %*% v)

ans.d <- diversitree:::expmv.expokit.dense(Q2, t, v)
expect_that(ans.e, equals(ans.d))

Q2.sparse <- diversitree:::expm.expokit.sparse.pars(Q2)
ans.s <- c(diversitree:::expmv.expokit.sparse(Q2, t, v, tol=1e-10))
expect_that(ans.s, equals(ans.d))

## Then, with a vector of times:
tt <- seq(0, t, length.out=11)[-1]

ans.ee <- sapply(tt, function(t) c(expm(Q2*t) %*% v))
ans.dd <- diversitree:::expmv.expokit.dense(Q2, tt, v)
ans.ss <- diversitree:::expmv.expokit.sparse(Q2.sparse, tt, v, tol=1e-10)

expect_that(ans.ee, equals(ans.dd))
expect_that(ans.ee, equals(ans.ss))
