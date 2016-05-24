
context("Test on ecd(2,-5)")

eps <- 0.001 # default tolerance of error for real number

gamma <- -5
d <- ecd(2,-5)
rnd <- rec(40000, quantilize(d))

test_that("test discriminant",{
  expect_true(discr(d) == -16*(-4*125+27*4))
})

test_that("test j-invariant",{
  expect_true(jinv(d) == 1728*20^3/discr(d))
})

test_that("test adj_gamma",{
    r2 <- ecd.adj2gamma(ecd.adj_gamma(gamma))
    expect_true(abs(r2-gamma) < eps)
})

test_that("test half of cdf",{
    c <- ecd.cdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test half of ccdf",{
    c <- ecd.ccdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test -Inf of cdf",{
    c <- ecd.cdf(d,-Inf)
    expect_true(abs(c-0) < eps)
})

test_that("test Inf of ccdf",{
    c <- ecd.ccdf(d,Inf)
    expect_true(abs(c-0) < eps)
})

test_that("test ccdf+cdf=1",{
    x <- 2
    c <- ecd.cdf(d, x)
    c0 <- ecd.ccdf(d, x)
    expect_true(abs(c+c0-1) < eps)
})

for(x in seq(-5,5,by=1)) {
    test_that(paste("test solve_sym at x=", x),{
        y1 <- solve(d,x)
        y2 <- solve_sym(d,x)
        expect_true(abs(y1-y2) < eps)
    })
}

test_that("test y_slope",{
    x <- 5
    dx <- 0.001
    dy <- y_slope(d,x)
    dy2 <- (solve(d,x+dx) - solve(d,x))/dx
    expect_true(abs(dy-dy2) < eps)
})

test_that("test ecd.pdf vs dec",{
    x <- 2
    p1 <- ecd.pdf(d,x)
    p2 <- dec(x, d)
    expect_true(abs(p1-p2) < eps)
})

test_that("test ecd.cdf vs pec",{
    x <- 2
    p1 <- ecd.cdf(d,x)
    p2 <- pec(x, d)
    expect_true(abs(p1-p2) < eps)
})

test_that("test mean of rec(40000)",{
    m <- mean(rnd)
    expect_true(abs(m-d@stats$mean) < 0.1)
})

test_that("test var of rec(40000)",{
    v <- var(rnd)
    expect_true(abs(v/d@stats$var-1) < 0.1)
})


