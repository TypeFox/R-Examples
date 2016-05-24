suppressMessages(library(Rmpfr))

context("Test on MPFR QAGI via J=0 Positive Dist")

# gamma=0, alpha>0
# x^2 = -y^3+alpha

eps <- 1e-3

y <- function(x, alpha) {
    xa <- alpha-x^2
    sign(xa)*abs(xa)^(1/3)
}
pdf <- function(x, alpha) exp(y(x, alpha)) 
pdf_inf <- function(t, alpha) pdf((1-t)/t, alpha)/t^2
pdf_minf <- function(t, alpha) pdf((t-1)/t, alpha)/t^2

tzero <- mpfr(.Machine$double.eps, 120L)

dh <- ecd(alpha=100)
dh2 <- ecd(alpha=100, sigma=ecd.mpfr(1), bare.bone=TRUE)

test_that("test const by integrate",{
    I <- integrate(pdf, 0, Inf, alpha=dh@alpha)
    expect_true(I$message=="OK" & abs(I$value*2/dh@const-1) < eps)
})

test_that("test const by mpfr_qagi",{
    tol <- 1e-4*dh@const
    I2 <- ecd.mpfr_qagi(dh2, pdf, 0, Inf, alpha=dh2@alpha, abs.tol=tol)
    expect_true(I2$message=="OK" & abs(I2$value*2/dh@const-1) < eps)
})

test_that("test cdf=1 by integrateR QAGI",{
    tol <- 1e-4*dh@const
    R2 <- integrateR(pdf_inf, tzero, 1, alpha=100, abs.tol=tol)
    expect_true(R2$message=="OK" & abs(R2$value*2/dh@const-1) < eps)
})

test_that("test eqv of integrate, integrateR QAGI at x=2",{
    tol <- 1e-4*dh@const
    # x has to be positive
    x <- 2
    t <- 1/(1+x)
    v1 <- integrate(pdf, x, Inf, alpha=100)$value
    v2 <- integrate(pdf_inf, 0, t, alpha=100)$value
    v3 <- integrateR(pdf_inf, tzero, t, alpha=100, abs.tol=tol)$value
    expect_true(abs(v2/v1-1) + abs(v3/v1-1) < eps)
})

test_that("test eqv of integrate, integrateR QAGI at x=-2",{
    tol <- 1e-4*dh@const
    # x has to be negative
    x <- -2
    t <- 1/(1-x)
    v1 <- integrate(pdf, -Inf, x, alpha=100)$value
    v2 <- integrate(pdf_minf, 0, t, alpha=100)$value
    v3 <- integrateR(pdf_minf, tzero, t, alpha=100, abs.tol=tol)$value
    expect_true(abs(v2/v1-1) + abs(v3/v1-1) < eps)
})

# --------------------------------------------

test_that("test eqv of const",{
    dh2 <- ecd(alpha=100, sigma=ecd.mpfr(1))
    expect_true(abs(dh2@const/dh@const-1) < 2*eps)
})

test_that("test eqv of integrating from 1 to Inf",{
    tol <- 1e-4*dh@const
    I1 <- integrate(pdf, 1, Inf, alpha=100)$value
    I2 <- ecd.mpfr_qagi(dh2, pdf, 1, Inf, alpha=100, abs.tol=tol)$value
    expect_true(abs(I2/I1-1) < eps)
})

test_that("test eqv of integrating from -Inf to -1",{
    tol <- 1e-4*dh@const
    I1 <- integrate(pdf, -Inf, -1, alpha=100)$value
    I2 <- ecd.mpfr_qagi(dh2, pdf, -Inf, -1, alpha=100, abs.tol=tol)$value
    expect_true(abs(I2/I1-1) < eps)
})

# -------------------------------------------------------------
# test moments

dh3 <- ecd(alpha=100, sigma=ecd.mpfr(1), with.stats=FALSE)

test_that("test eqv of integrating 2nd moment from 2 to Inf",{
    C <- ecd.mp2f(dh3@const)
    tol <- 1e-4 #*dh3@const
    mnt2 <- function(x) pdf(x,100)/C*x^2
    I1 <- integrate(mnt2, 2, Inf)$value
    I2 <- ecd.mpfr_qagi(dh3, mnt2, 2, Inf)$value
    expect_true(abs(I2/I1-1) < eps)
})

test_that("test eqv of integrating 4th moment from 2 to Inf",{
    C <- ecd.mp2f(dh3@const)
    tol <- 1e-4 #*dh3@const
    mnt4 <- function(x) pdf(x,100)/C*x^4
    I1 <- integrate(mnt4, 2, Inf)$value
    I2 <- ecd.mpfr_qagi(dh3, mnt4, 2, Inf)$value
    expect_true(abs(I2/I1-1) < eps)
})