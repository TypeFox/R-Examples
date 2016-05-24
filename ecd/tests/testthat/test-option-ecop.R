
suppressMessages(library(fOptions))

context("Test on Option Pricing")

eps <- 0.001 # default tolerance of error for real number

# two SPX options from 2015-03-19, a day before triple witching
C <- c(0.05, 1.8, 50, 140.2)
K <- c(2200, 2100, 2040, 1950)
S <- 2089.27
T <- 1/365
y <- 0.019
vol0 <- c(0.3710537, 0.128886, 0.294296, 0.7237632) 

# ------------------------------------------------------
test_that("test the BS implied volatility",{
    vol <- ecop.bs_implied_volatility(C, K, S, T, div_yield=y)
    expect_true(sum(abs(vol/vol0-1)) < eps)
})

test_that("test the BS implied volatility using MPFR",{
    vol <- ecop.bs_implied_volatility(C*ecd.mp1, K, S, T, div_yield=y)
    expect_true(sum(abs(vol/vol0-1)) < eps)
})

# ------------------------------------------------------
ki <- seq(-10, 10)

test_that("test k1 formula used in OGF star",{
    sigma = 0.001
    mu = 0.002

    k <- ki*sigma+mu
    
    sigma1 = 0.003
    mu_D1 = -sigma1^2/4
    
    ki1a <- sigma/sigma1*ki + (mu/sigma1+sigma1/4)
    ki1b <- (k-mu_D1)/sigma1
    
    expect_true(max(abs(ki1a-ki1b)) < eps)
})

test_that("test exact formula for lambda=1 OGF star",{
    ld1 <- ecld(lambda=1)

    L1 <- ecld.ogf_star(ld1, ki)
    L2 <- 1/2/sqrt(pi)*(exp(-ki^2) - sqrt(pi)*abs(ki)*Rmpfr::erfc(abs(ki)))
    
    expect_true(max(abs(L1/L2-1)) < eps)
})

# ------------------------------------------------------
test_that("test asymptotic formula for lambda=1 OGF star",{
    ki <- c(5, 8, 10)
    ld1 <- ecld(lambda=1)
    
    L1 <- ecld.ogf_star(ld1, ki)
    L2 <- ecld.ogf_star_exp(ld1, ki, order=4)
    
    expect_true(max(abs(L1/L2-1)) < eps)
})
