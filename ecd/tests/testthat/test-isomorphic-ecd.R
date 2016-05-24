
context("Test on Isomorphic Mapping")

eps <- 0.001 # default tolerance of error for real number

theta <- 0.25*pi # 45 degree

# on unit circle
d0 <- ecd.polar(R=1,theta=theta)
a0 <- d0@alpha
r0 <- d0@gamma

x <- 0.5
y0 <- solve(d0, x)

# iso mapping
lambda <- 2
s1 <- d0@sigma*lambda^3

a1 <- a0*lambda^6
r1 <- r0*lambda^4
d1 <- ecd(a1,r1)
d1s <- ecd(a1,r1,1/s1)

# sigma adjusted
a2 <- a0*s1^2
r2 <- r0*s1^(4/3)

test_that("test a1 = a2",{
    expect_true(abs(a1-a2) < eps)
})

test_that("test r1 = r2",{
    expect_true(abs(r1-r2) < eps)
})

test_that("test j-inv of d0 and d1",{
    expect_true(abs(jinv(d0)-jinv(d1)) < eps)
})

test_that("test s1",{
    expect_true(abs(jinv(d0)-jinv(d1)) < eps)
})

test_that("test y(x) scaled by lambda",{
    y1 <- solve(d0, x)*lambda^2
    y2 <- solve(d1, x*lambda^3)
    y3 <- solve(d1s, x)
    y4 <- solve(d1, x*s1)
    delta <- abs(y2-y1) + abs(y3-y1) + abs(y4-y1)
    expect_true(delta < eps)
})

test_that("test y(x) scaled by sigma and R",{
    y1 <- solve(d1, x)
    y2 <- solve(d1s, x/s1)
    y3 <- solve(d0, x/s1)*s1^(2/3)
    y4 <- solve(d0, x/sqrt(d1@R))*d1@R^(1/3)
    delta <- abs(y2-y1) + abs(y3-y1) + abs(y4-y1)
    expect_true(delta < eps)
})

# special case for y(0), scales at R^(1/3)
test_that("test y(0) scaling",{
    R <- 30^3
    t <- 45/180*pi
    y1 <- solve(ecd.polar(R=1, theta=t, bare.bone=TRUE), 0)
    y2 <- solve(ecd.polar(R=R, theta=t, bare.bone=TRUE), 0)
    expect_true(abs(y2/y1-R^(1/3)) < eps)
})

# y0 isomorphism
for (R in c(1, 8, 27)) {
    for (degree in c(0, 35, 90, 180, 210, 270, 310, 315)) {
        test_that(paste("test y(0) isomorphism at R=",R,"deg=",degree),{
            theta <- degree/180*pi
            d <- ecd.polar(R,theta, bare.bone=TRUE)
            y0 <- solve(d,0)
            y1 <- ecd.y0_isomorphic(theta, d@R)
            y2 <- ecd.y0_isomorphic(object=d)
            expect_true(abs(y1-y0) + abs(y2-y0) < eps)
        })
    }
}

# acosh(1/tan) causes NaN between 180 and 225 sporadically
test_that("test NaN of y(0) isomorphism",{
    degree <- seq(0, 315, by=0.01)
    y0 <- ecd.y0_isomorphic(degree/180*pi)
    expect_true(all(!is.na(y0)))
})

# -------------------------------------------------
# j-inv in polar form

test_that("test j-inv, r_adj>0",{
	t <- 30/180*pi
    j1 <- jinv(ecd.polar(R=5, theta=t))
    j2 <- 1728*sin(t)^2
    expect_true(abs(j2-j1) < eps)
})

test_that("test j-inv, r_adj<0",{
    t <- (180+30)/180*pi
    j1 <- jinv(ecd.polar(R=5, theta=t))
    j2 <- 1728/(1-tan(t)^(-2))
    expect_true(abs(j2-j1) < eps)
})


