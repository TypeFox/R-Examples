
context("Test on ecd(-1,-1)")

dx <- 0.01
x <- seq(-100, 100, by=dx)
eps <- 0.001 # default tolerance of error for real number

d <- ecd(-1,-1, with.quantile=TRUE)

d1C <- ecd(-1,ecd.adj2gamma(-1))
d1R <- ecd.polar(R=sqrt(2), theta=225/180*pi, with.quantile=TRUE)

test_that("test known const",{
    c <- 1.6666
    expect_true(abs(log(const(d)/c)) < eps)
})

test_that("test theta on ecd(-1,ecd.adj2gamma(-1))",{
    expect_true(abs(d1C@theta/pi*180-225) < eps)
})

test_that("test theta by ecd.polar",{
    expect_true(abs(d1R@theta/pi*180-225) < eps)
})

test_that("test R by ecd.polar",{
    expect_true(abs(d1R@R-sqrt(2)) < eps)
})

test_that("test jinv of ecd.polar",{
    j1 <- jinv(d)/1728
    j2 <- (1-1/tan(d@theta)^2)^(-1)
    expect_true(abs(j1-j2) < eps)
})

test_that("test jinv of ecd(-1,1D)",{
    d2 <- ecd(-1,ecd.adj2gamma(1))
    j1 <- jinv(d2)/1728
    j2 <- sin(135/180*pi)^2
    expect_true(abs(j1-j2) < eps)
})

test_that("test equality of kurtosis between d1R and d1C",{
    expect_true(abs(d1R@stats$kurtosis - d1C@stats$kurtosis) < eps)
})

test_that("test half of cdf",{
    c <- ecd.cdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test half of cdf for d1R",{
    c <- ecd.cdf(d1R,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test discriminant",{
  expect_true(discr(d) == -16*23)
})

test_that("test j-invariant",{
  expect_true(jinv(d) == 1728*4^3/discr(d))
})

test_that("test y_slope",{
    d <- ecd(-1, -1, 0.5, 0.1, 0.1)
    x <- 2
    dx <- 0.001
    slope <- (solve(d,x+dx)-solve(d,x))/dx
    slope2 <- y_slope(d,x)
    expect_true(abs(slope-slope2) < eps)
})
test_that("test analytic y_slope",{
    b <- 0.1
    d <- ecd(-1, -1, 1, b)
    x <- 2
    y <- solve(d,x)
    slope <- -(2*x+b*y)/(3*y^2+(-1)+b*x)
    slope2 <- y_slope(d,x)
    expect_true(abs(slope-slope2) < eps)
})


test_that("test ellipticity",{
    xe <- ellipticity(d)$avg
    xe2 <- 2.2767
    expect_true(abs(xe-xe2) < eps)
})

test_that("test ellipticity on ecd(-27,0)",{
    d <- ecd(-27,0)
    xe <- ellipticity(d)$avg
    xe2 <- 9
    expect_true(abs(xe/xe2-1) < eps)
})

test_that("test ellipticity on ecd(4,0)",{
    d <- ecd(4,0)
    xe <- ellipticity(d)$avg
    xe2 <- 2
    expect_true(abs(xe/xe2-1) < eps)
})


for(x in seq(-5,5,by=1)) {
    test_that(paste("test solve_sym at x=", x),{
        y1 <- solve(d,x)
        y2 <- solve_sym(d,x)
        expect_true(abs(y1-y2) < eps)
    })
}

test_that("test asymp moment(0) against cdf",{
    x <- 10
    c1 <- 1-moment(d, 0, asymp.lower=-x, asymp.upper=x)
    c2 <- ecd.cdf(d,-x) + ecd.ccdf(d,x)
    expect_true(abs(c1-c2) < eps)
})

test_that("test asymp kurtosis",{
    k1 <- ecd.asymp_kurtosis(d,0)
    k2 <- ecd.kurtosis(d)
    expect_true(abs(k1-k2) < eps)
})
