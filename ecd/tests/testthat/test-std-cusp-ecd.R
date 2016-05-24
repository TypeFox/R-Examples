
context("Test on Stanrdard Cusp Distribution")

eps <- 0.001 # default tolerance of error for real number

d <- ecd()
stats <- d@stats
c0 <- 1.5*sqrt(pi)

test_that("test the model short name",{
    expect_true(d@model[3] == "symm std csp ecd")
})

test_that("test known const",{
  expect_true(abs(const(d)-c0) < eps)
})

test_that("test R=0 and theta=0",{
  expect_true(d@R^2+d@theta^2 == 0)
})

test_that("test half of cdf",{
    c <- ecd.cdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test half of ccdf",{
    c <- ecd.ccdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test mu2",{
  m2c <- 13.125
  expect_true(abs(log(stats$m2/m2c)) < eps)
})

test_that("test odd std moments being zero",{
    n <- seq(1,11,by=2)
    m <- ecd.cusp_std_moment(n)
    expect_true(sum(m)==0)
})
test_that("test first six even std moments",{
    n <- seq(0,6,by=2)
    m <- ecd.cusp_std_moment(n)
    m2 <- c(1, 13.125, 2111.484, 1278767.725)
    expect_true(sum(abs(log(m/m2))) < eps)
})

test_that("test kurtosis",{
  kt <- stats$kurtosis
  kt0 <- 429/35
  expect_true(abs(log(kt/kt0)) < eps)
})

test_that("test effect of sigma on const",{
    s <- 0.5
    d <- ecd(0,0,s)
    p <- ecd.pdf(d,1)
    p1 <- exp(-(1/s)^(2/3))/s/c0
    expect_true(abs(log(p/p1)) < eps)
})

test_that("test with.stats = FALSE",{
    d2 <- ecd(with.stats=FALSE)
    expect_true(length(d2@stats)==0)
})

test_that("test effect of mu on mean",{
    mu <- 1
    d <- ecd(0,0,1,0,mu)
    m1 <- d@stats$mean
    expect_true(abs(m1-mu) < eps)
})

test_that("test sum on segments of cdf",{
    x0 <- c(-Inf, -10, 0, 10)
    x <- ecd.lag(x0,-1)
    x[length(x0)] <- Inf
    c <- sum(ecd.cdf(d, x, from.x = x0))
    expect_true(abs(c-1) < eps)
})

test_that("test piece.wise mode of cdf",{
    x <- seq(-100,10,by=0.5)
    c <- ecd.cdf(d, x, piece.wise=TRUE)
    cp <- c[length(x)]
    c0 <- ecd.cdf(d, max(x))
    expect_true(abs(cp-c0) < eps)
})

test_that("test 2 vectors of cdf",{
    x <- c(0,Inf)
    c <- sum(ecd.cdf(d, x))
    expect_true(abs(c-1.5) < eps)
})

for (x in seq(-20, -5, by=5)) {
	test_that(paste("test cdf vs pgamma at x=",x),{
		c1 <- ecd.cdf(d, x)
		t <- abs(x)^(2/3)
		c2 <- (1-pgamma(t,shape=3/2))/2 
		expect_true(abs(c1-c2) < eps)
	})
}

for (x in seq(1, 13, by=3)) {
	test_that(paste("test cdf vs erfc at x=",x),{
		c1 <- ecd.ccdf(d, x)
        erfc <- function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
		c2 <- 1/sqrt(pi)*x^(1/3)*exp(-x^(2/3)) + 1/2*erfc(x^(1/3))
		expect_true(abs(c2/c1-1) < eps)
	})
}

test_that("test asymptotic expansion of upper incomplete gamma",{
    x <- -20
    c1 <- ecd.cdf(d, x)
    
    t <- abs(x)^(2/3)
	k <- seq(0, 10)
	fr <- function(k) { gamma(3/2)/gamma(3/2-k)/t^k } 
	S <- sapply(k, fr)
	c2 <- t^(1/2)*exp(-t)*sum(S)/sqrt(pi)

    expect_true(abs(c2/c1-1) < eps)
})

test_that("test piece.wise mode of ccdf",{
    x <- seq(-10,10,by=1)
    c <- ecd.ccdf(d, x, piece.wise=TRUE)
    c0 <- ecd.cdf(d, x, piece.wise=FALSE)
    expect_true(abs(sum(c-c0)) < eps)
})

test_that("test ccdf+cdf=1",{
    x <- 5
    c <- ecd.cdf(d, x)
    c0 <- ecd.ccdf(d, x)
    expect_true(abs(c+c0-1) < eps)
})

test_that("test imgf",{
	x <- 0.03
	s <- 0.05
	d1 <- ecd(0,0, s)
	d2 <- ecd(0,0)
	I1 <- ecd.imgf(d1,x)
	I2 <- ecd.imgf(d2, x/s, t=s)
    expect_true(abs(I1-I2) < eps)
})

test_that("test option drift",{
	s <- 0.05
	d1 <- ecd(0,0, s)
	d2 <- ecd(0,0)
	I1 <- ecd.imgf(d1)
	I2 <- ecd.imgf(d2, t=s)
    expect_true(abs(I1-I2) < eps)
})

test_that("test one absolute value in std cusp analytic mgf",{
    M1 <- 1.010733
    M2 <- ecd.cusp_std_mgf(1, sigma=0.04)
    expect_true(abs(M1-M2) < eps)
})

test_that("test one absolute value in std cusp imgf",{
    M1 <- 1.010733
    M2 <- ecd.imgf(ecd(sigma=0.04))
    expect_true(abs(M1-M2) < eps)
})

test_that("test mu in std cusp analytic mgf",{
    t <- 1.2
    M0 <- ecd.cusp_std_mgf(t, sigma=0.04)
    M1 <- ecd.cusp_std_mgf(t, mu=1, sigma=0.04)
    expect_true(abs(M1/M0-exp(t)) < eps)
})

for(s in seq(0.01, 0.04, by=0.01)) {
    mu1 <- -log(ecd.imgf(ecd(), -Inf, t=s))
    
    test_that(paste("test mu_D scaling at sigma", s),{
        mu2 <- -log(ecd.imgf(ecd(sigma=s), -Inf))
        expect_true(abs(exp(mu1-mu2)-1) < eps)
    })
    test_that(paste("test mu_D with std mgf at sigma", s),{
        mu3 <- -log(ecd.cusp_std_mgf(s))
        expect_true(abs(exp(mu1-mu3)-1) < eps)
    })

    # non-zero mu is required to test edge cases
    M1 <- try(ecd.imgf(ecd(sigma=s, mu=1), -Inf, unit.sigma=FALSE))

    test_that(paste("test unit.sigma of imgf at sigma", s),{
        M2 <- ecd.imgf(ecd(sigma=s, mu=1), -Inf, unit.sigma=TRUE)
        expect_true(abs(M2/M1-1) < eps)
    })
    test_that(paste("test risk neutral isomorphism at sigma", s),{
        M3 <- ecd.imgf(ecd(mu=1/s), -Inf, t=s)
        expect_true(abs(M3/M1-1) < eps)
    })   
    test_that(paste("test imgf vs std cusp analytic form at sigma", s),{
        M4 <- ecd.cusp_std_mgf(1, mu=1, sigma=s)
        expect_true(abs(M4/M1-1) < eps)
    })
}



