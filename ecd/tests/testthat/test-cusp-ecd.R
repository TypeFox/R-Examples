
context("Test on General Cusp Distribution")

eps <- 0.001 # default tolerance of error for real number

# test ecd.devel
test_that("test ecd.devel is set",{
    dev <- ecd.devel()
    expect_true(dev == TRUE || dev == FALSE)
})

#
d <- ecd.cusp(alpha=1)
d2 <- ecd.cusp(gamma=-1)

test_that("test the model short name, type 1",{
    expect_true(d@model[3] == "symm csp ecd")
})

test_that("test the model short name, type 2",{
    expect_true(d2@model[3] == "symm csp ecd")
})

test_that("test d: theta is 315 degree",{
    expect_true(abs(d@theta/pi*180-315) < eps)
})

test_that("test d2: theta is 315 degree",{
    expect_true(abs(d2@theta/pi*180-315) < eps)
})

test_that("test half of cdf",{
    c <- ecd.cdf(d,0)
    expect_true(abs(c-0.5) < eps)
})

test_that("test effect of mu on mean",{
    mu <- 1
    d <- ecd.cusp(alpha=1, mu=mu)
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

# testing the stability of cubic solution 
for(alpha in c(1, 5, 10, seq(20, 100,by=20), 10^c(3,4,5,6,7))) {
    test_that(paste("test solve at x=0 for alpha",alpha),{
        d <- ecd.cusp(alpha=alpha)
        y1 <- solve(d, 0)
        y2 <- min(Re(ecd.cubic(d,0)[[1]]))
        expect_true(abs(y1-y2) < eps)
    })
}

test_that("test y(0) of cusp",{
    y0 <- solve(d2,0)
    y1 <- -(d2@alpha/2)^(1/3)
    expect_true(abs(y0-y1) < eps)
})

# ---------------------------------------------
# asym cusp

test_that("test trig solution of asym cusp for beta=0.5",{
    d <- ecd(beta=0.5)
    x0 <- -(4*d@beta^3)/27 # this is where the scenario changes
    x <- seq(-abs(x0)*10, abs(x0)*10, length.out=100)
    y1 <- solve(d,x)
    y2 <- ecd.solve_cusp_asym(x,d@beta)

    expect_true(sum(abs(y1-y2)) < eps)
})

test_that("test trig solution of asym cusp for beta=-0.6",{
    d <- ecd(beta=-0.6)
    x0 <- -(4*d@beta^3)/27 # this is where the scenario changes
    x <- seq(-abs(x0)*10, abs(x0)*10, length.out=100)
    y1 <- solve(d,x)
    y2 <- ecd.solve_cusp_asym(x,d@beta)
    
    expect_true(sum(abs(y1-y2)) < eps)
})

beta <- seq(-0.8, 0.8, by=0.05)
beta <- beta[beta != 0]
dasym <- parallel::mclapply(beta, function(b) ecd(beta=b))    

test_that("test lm of m1 and skewness of asym cusp",{
    m1 <- sapply(dasym, function(d) d@stats$m1)
    S <- sapply(dasym, function(d) d@stats$skewness)

    # fit1 <- lm(m1 ~ 0+beta)
    # fit2 <- lm(S ~ 0+beta)
    # fit3 <- lm(m1S ~ poly(S,2,raw=TRUE))
    # plot (S, m1S)
    # lines (S, predict(fit3,S=data.frame(S)))
    m1p <- 1.2353 * beta
    Sp <- 1.2113 * beta
    m1Sp <- sign(S)* (0.00178 + 0.9875 *abs(S) + 0.0405 * abs(S)^2)   
    # --------------------------
    err1 <- max(abs(m1p/m1-1))
    err2 <- max(abs(Sp/S-1))
    err3 <- max(abs(m1Sp/m1-1))
    
    err <- max(c(err1, err2, err3))
    expect_true(err < 0.03)
})

# test y+ and y-
fasym_C <- function(x,beta) {
    y1 <- ecd.solve_cusp_asym(x,beta)
    y2 <- ecd.solve_cusp_asym(-x,beta)
    exp(y1)+exp(y2)
}
fasym_mn <- function(x,beta,n) {
    y1 <- ecd.solve_cusp_asym(x,beta)
    y2 <- ecd.solve_cusp_asym(-x,beta)
    x^n*(exp(y1)+(-1)^n*exp(y2))
}

for (beta in c(0.6, 0.5, 0.4)) {
    d <- ecd(beta=beta)
    C <- integrate(fasym_C, 0, Inf, beta=beta)$value

    test_that(paste("test C of asym cusp by y+/- at beta=", beta),{
        expect_true(abs(C/d@const-1) < eps)
    })
    
    m1 <- integrate(function(x) fasym_mn(x,beta,1), 0, Inf)$value/C
    m2 <- integrate(function(x) fasym_mn(x,beta,2), 0, Inf)$value/C
    m3 <- integrate(function(x) fasym_mn(x,beta,3), 0, Inf)$value/C

    test_that(paste("test m1, m2, m3 of asym cusp by y+/- at beta=", beta),{
        e1 <- abs(m1/d@stats$m1-1) 
        e2 <- abs(m2/d@stats$m2-1)
        e3 <- abs(m3/d@stats$m3-1)
        expect_true(e1+e2+e3 < eps)
    })
    
    test_that(paste("test B+ - B- law at beta=", beta),{
        a <- ecd.cdf(d, 0)*d@const
        b <- ecd.ccdf(d, 0)*d@const
        e1 <- abs(b-a - d@beta)
        expect_true(e1 < eps)
    })

}
