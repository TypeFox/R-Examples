context("Tests of Beta Binomial and Chance-Corrected Beta-Binomial log-likelihood functions")


x <- c(2, 2, 3, 4, 5, 5, 5, 5, 6, 6, 6, 7,
       7, 7, 7, 7, 7, 8, 8, 8, 9, 11, 12)
X2 <- cbind(x, 12)

summary(fit <- betabin(X2, corrected=FALSE, method="duotrio"))
summary(fita <- betabin(X2, corrected=FALSE, method="triangle"))

summary(fit2 <- betabin(X2, corrected=TRUE, method="duotrio"))
summary(fit2a <- betabin(X2, corrected=TRUE, method="triangle"))

#################################
## Some functions:
mg2ab <- function(par) {
    par <- unname(par)
    a <- par[1] * (1 - par[2]) / par[2]
    b <- (1 - par[2]) * (1 - par[1]) / par[2]
    c(alpha=a, beta=b)
}

BB <- function(X, a, b) {
    fun <- function(x, n) {
        lchoose(n, x) - lbeta(a, b) + lbeta(a + x, b + n - x)
    }
    sum(apply(X, 1, function(xx) fun(xx[1], xx[2])))
}

CBB <- function(X, a, b, pg) {
    fun <- function(x, n) {
        SUM <- 0
        for(i in 0:x) {
            SUM <- SUM + choose(x, i) *
                (1 - pg)^(n - x + i) * pg^(x - i) *
                    beta(a + i, n - x + b)
        }
        lchoose(n, x) - lbeta(a, b) + log(SUM)
    }
    sum(apply(X, 1, function(xx) fun(xx[1], xx[2])))
}
#################################

## Evaluate likelihood of beta-binomial model
## 1) at optimum for pg=1/2:
## 2) at optimum for pg=1/3:
## 3) at (mu, gamma) = (.5, .5) for pg=1/2
## 4) at (mu, gamma) = (.5, .5) for pg=1/3

test_that("Beta-binomial log-likelihood evaluation is correct", {

    ## 1)
    (AB <- mg2ab(coef(fit)))
    (Y <- BB(X2, AB[1], AB[2]))
    (Z <- logLik(fit))
    expect_equal(Y, as.vector(Z))

    ## 2)
    (AB <- mg2ab(coef(fita)))
    (Y <- BB(X2, AB[1], AB[2]))
    (Z <- logLik(fita))
    expect_equal(Y, as.vector(Z))

    ## 3)
    rho <- update(fit, doFit=FALSE)
    rho$par
    (Y <- rho$Factor - sensR:::setParBB(rho))
    (AB <- mg2ab(rho$par))
    (Z <- BB(X2, AB[1], AB[2]))
    expect_equal(Y, Z)

    ## 4)
    rhoa <- update(fita, doFit=FALSE)
    rhoa$par
    (Y <- rhoa$Factor - sensR:::setParBB(rhoa))
    (AB <- mg2ab(rhoa$par))
    (Z <- BB(X2, AB[1], AB[2]))
    expect_equal(Y, Z)
})
#################################

## Evaluate likelihood of *chance-corrected* beta-binomial model
## 1) at optimum for pg=1/2:
## 2) at optimum for pg=1/3:
## 3) at (mu, gamma) = (.5, .5) for pg=1/2
## 4) at (mu, gamma) = (.5, .5) for pg=1/3

test_that("*Chance-Corrected* Beta-binomial log-lik eval is correct", {
    ## 1)
    AB <- mg2ab(coef(fit2))
    (Y <- CBB(X2, AB[1], AB[2], pg=1/2))
    (Z <- logLik(fit2))
    expect_equal(Y, as.vector(Z))

    ## 2)
    AB <- mg2ab(coef(fit2a))
    (Y <- CBB(X2, AB[1], AB[2], pg=1/3))
    (Z <- logLik(fit2a))
    expect_equal(Y, as.vector(Z))

    ## 3)
    rho2 <- update(fit2, doFit=FALSE)
    rho2$par
    (Y <- rho2$Factor - sensR:::setParBB(rho2))
    (AB <- mg2ab(rho2$par))
    (Z <- CBB(X2, AB[1], AB[2], pg=1/2))
    expect_equal(Y, Z)

    ## 4)
    rho2a <- update(fit2a, doFit=FALSE)
    rho2$par
    (Y <- rho2a$Factor - sensR:::setParBB(rho2a))
    (AB <- mg2ab(rho2a$par))
    (Z <- CBB(X2, AB[1], AB[2], pg=1/3))
    expect_equal(Y, Z)
})
