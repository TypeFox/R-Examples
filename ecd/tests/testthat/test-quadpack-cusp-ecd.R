suppressMessages(library(Rmpfr))

context("Test on MPFR QAGI via Cusp Dist")

eps <- 1e-4

c0 <- 2/3/sqrt(pi)
pdf <- function(x) exp(-(x^2)^(1/3)) * c0
pdf_inf <- function(t) pdf((1-t)/t)/t^2
pdf_minf <- function(t) pdf((t-1)/t)/t^2

tzero <- mpfr(.Machine$double.eps, 120L)

test_that("test cdf=1 by integrate",{
    I <- integrate(pdf, 0, Inf)
    expect_true(I$message=="OK" & abs(I$value*2-1) < eps)
})

test_that("test cdf=1 by integrate QAGI",{
    I2 <- integrate(pdf_inf, 0, 1)
    expect_true(I2$message=="OK" & abs(I2$value*2-1) < eps)
})

test_that("test cdf=1 by integrateR QAGI",{
    R2 <- integrateR(pdf_inf, tzero, 1)
    expect_true(R2$message=="OK" & abs(R2$value*2-1) < eps)
})

test_that("test eqv of integrate, integrateR QAGI at x=2",{
    # x has to be positive
    x <- 2
    t <- 1/(1+x)
    v1 <- integrate(pdf, x, Inf)$value
    v2 <- integrate(pdf_inf, 0, t)$value
    v3 <- integrateR(pdf_inf, tzero, t)$value
    expect_true(abs(v2/v1-1) + abs(v3/v1-1) < eps)
})

test_that("test eqv of integrate, integrateR QAGI at x=-2",{
    # x has to be negative
    x <- -2
    t <- 1/(1-x)
    v1 <- integrate(pdf, -Inf, x)$value
    v2 <- integrate(pdf_minf, 0, t)$value
    v3 <- integrateR(pdf_minf, tzero, t)$value
    expect_true(abs(v2/v1-1) + abs(v3/v1-1) < eps)
})

# --------------------------------------------
d0 <- ecd()
d1 <- ecd(sigma=ecd.mpfr(1), with.stats=FALSE)

test_that("test eqv of const",{
    expect_true(abs(d1@const/d0@const-1) + abs(d1@const*c0-1) < 2*eps)
})

test_that("test eqv of integrating from 1 to Inf",{
    I1 <- integrate(pdf, 1, Inf)$value
    I2 <- ecd.mpfr_qagi(d1, pdf, 1, Inf)$value
    expect_true(abs(I2/I1-1) < eps)
})

test_that("test eqv of integrating from -Inf to -1",{
    I1 <- integrate(pdf, -Inf, -1)$value
    I2 <- ecd.mpfr_qagi(d1, pdf, -Inf, -1)$value
    expect_true(abs(I2/I1-1) < eps)
})

# --------------------------------------------
test_that("test mpfr_qagi by eqv of integrating from 1 to Inf of small sigma",{
    pdf <- function(x,s) exp(-((x/s)^2)^(1/3)) 

    s3 <- 0.01
    d3 <- ecd(sigma=s3)
    I3 <- integrate(pdf, s3, Inf, s=s3)$value/s3

    s4 <- 0.00001
    d4 <- ecd(sigma=ecd.mpfr(s4), bare.bone=TRUE)
    I4 <- ecd.mpfr_qagi(d4, pdf, s4, Inf, s=s4)$value/s4
    expect_true(abs(I3/I4-1) < eps)
})

