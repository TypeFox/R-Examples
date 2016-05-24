
context("Test on Various Taylor Expansions")

eps <- 0.001 # default tolerance of error for real number

# ---------------------------------------------
# cusp y(z)

dc <- ecd.cusp(alpha=0.5, sigma=2)
x <- 0.25
z <- sqrt(x^2/(dc@alpha*dc@sigma^2))
y_x <- solve(dc,x)

test_that("test y_x in trig form",{
    y_x_trig <- solve_trig(dc,x)
    expect_true(abs(y_x_trig/y_x-1) < eps)
})

test_that("test y_z in trig form",{
    y_z <- -(4*dc@alpha)^(1/3)*cos(1/3*acos(z^2-1))
    expect_true(abs(y_z/y_x-1) < eps)
})

test_that("test y_z in series form",{
    y_z_pow <- -(dc@alpha/2)^(1/3)*(1 + sqrt(2/3)*z -z^2/9 
                                    + 5*z^3/2^(3/2)/3^(7/2)
                                    + 4*z^4/243)
    expect_true(abs(y_z_pow/y_x-1) < eps)
})

a <- 10^5
s <- 1
dc2 <- ecd.cusp(alpha=a, sigma=ecd.mpfr(s))

test_that("test asymptotics to Laplace distribution",{
    x <- 0.5
    y_x <- solve(dc2,x)
    y_x2 <- -(a/2)^(1/3) - (2/a)^(1/6)*x/(sqrt(3)*s)
    expect_true(abs(y_x2/y_x-1) < eps)
})

test_that("test var of asymptotics for cusp",{
    v0 <- ecd.sd(dc2)^2/s^2
    v1 <- (10 + 3*2^(2/3)*a^(1/3) + 13/(3*2^(2/3)*a^(1/3)) 
           -17*2^(2/3)/(9*a^(2/3)) + 445/(108*a))
    expect_true(abs(v1/v0-1) < eps)
})

test_that("test kurtosis of asymptotics for cusp",{
    k0 <- ecd.kurtosis(dc2)
    k1 <- (6 + 2^(10/3)/a^(1/3) - 2^(8/3)/a^(2/3) - 44/(3*a)
           + 19*2^(10/3)/(3*a^(4/3)) - 257*2^(2/3)/(3*a^(5/3))
           - 15140285/(27*a^2))
    expect_true(abs(k1/k0-1) < eps)
})

# ---------------------------------------------
# j=0

a <- 10^3
s <- 1
d1j0 <- ecd(alpha=a, sigma=ecd.mpfr(s))
d2j0 <- ecd(alpha=-a, sigma=ecd.mpfr(s))

test_that("test var of asymptotics for j=0",{
    f_var <- function(a) {
        a23 <- (a^2)^(1/3)
        a13 <- sign(a)*abs(a)^(1/3)
        ( 63/8 +3/2*a23 -9/2*a13 -39/8/a13 -117/32/a23 -3/4/a )
    }
    
    v1 <- ecd.sd(d1j0)^2/s^2
    v2 <- ecd.sd(d2j0)^2/s^2
    
    expect_true(abs(f_var(a)/v1-1) + abs(f_var(-a)/v2-1) < eps)
})

test_that("test kurtosis of asymptotics for j=0",{
    f_kurt <- function(a) {
        a23 <- (a^2)^(1/3)
        a13 <- sign(a)*abs(a)^(1/3)
        ( 3 -6/a13 +15/a23 -345/4/a23^2 +75/8/a13^5 )
    }
    
    v1 <- ecd.kurtosis(d1j0)
    v2 <- ecd.kurtosis(d2j0)
    
    expect_true(abs(f_kurt(a)/v1-1) + abs(f_kurt(-a)/v2-1) < eps)
})

# ---------------------------------------------
# j=1728

r <- 100
s <- 1
d1r <- ecd(gamma=r, sigma=ecd.mpfr(s))
d2r <- ecd(gamma=-r, sigma=ecd.mpfr(s))

test_that("test var of asymptotics for j=1728",{
    f_var <- function(r) {
        ( 45/8 +r/2 +135/4/r -13365/32/r^2 +103275/8/r^3 -403920675/512/r^4 )
    }
    
    v1 <- ecd.sd(d1r)^2/s^2
    # v2 <- ecd.sd(d2r)^2/s^2 # can't do this, there is no analytic continuity here!
    
    expect_true(abs(f_var(r)/v1-1) < eps)
})

test_that("test kurtosis of asymptotics for j=1728",{
    f_kurt <- function(r) {
        ( 3 +45/r +4725/4/r^2 -150255/16/r^3 -41802075/64/r^4 )
    }

    v1 <- ecd.kurtosis(d1r)
    # v2 <- ecd.kurtosis(d2r)

    expect_true(abs(f_kurt(r)/v1-1) < eps)
})

# ---------------------------------------------
# asymptotic variance

Ra <- 10000
ta <- pi/2
da <- ecd.polar(Ra, ta)

var_N <- function(d) {
    t <- d@theta
    T <- 1/3*log(tan(t/2))
    A <- 2^(5/3)/3*sin(t)^(1/3)*cosh(T)
    d@R^(2/3)*d@sigma^2/A
}

test_that("test asymp var",{
    s1 <- var_N(da)^(1/2)
    s2 <- ecd.sd(da)
    expect_true(abs(s2/s1-1) < 0.01)
})

test_that("test C estimate via asymp var",{
    c1 <- da@const
    c2 <- sqrt(2*pi*var_N(da))*exp(solve(da,0))
    expect_true(abs(c2/c1-1) < 0.01)
})

# test several representative angles
degree <- c(1, 45, 90, 135, 179) # singular on j=0
theta <- degree/180*pi
alpha <- Ra*cos(theta)
gamma <- ecd.adj2gamma(Ra*sin(theta))
y0 <- ecd.y0_isomorphic(theta, Ra)
var_N_cart <- (3*y0^2+gamma)/2

test_that("test cartesian expansion of y(0)",{
    y0_cart <- function(alpha, gamma) {
        d <- ecd(alpha, gamma, bare.bone=TRUE)
        phi <- (alpha+sqrt(-discr(d)/432))^(1/3)
        phi/2^(1/3) - 2^(1/3)*gamma/3/phi
    }
    y0c <- mapply(y0_cart, alpha, gamma, SIMPLIFY=TRUE)
    expect_true(max(abs(y0c-y0)) < eps)
})

test_that("test cartesian expansion of var_N with polar expansion",{
    d <- mapply(function(a,r) ecd(a,r,bare.bone=TRUE), alpha, gamma)
    var_N_polar <- sapply(d, function(d) var_N(d))
    expect_true(max(abs(var_N_cart/var_N_polar-1)) < eps)
})

test_that("test heuristic expansion on positive j=1728",{
    v1 <- var_N_cart[3]
    v2 <- gamma[3]/2
    expect_true(abs(v1/v2-1) < eps)
})

test_that("test heuristic expansion on j=0",{
    i <- c(1,5)
    y0j0 <- y0[i]
    y0err <- abs(alpha[i])^(1/3)/abs(y0j0)-1
    
    varj0 <- var_N_cart[i]
    varerr <- abs(alpha[i])^(2/3)*3/2/varj0-1
    
    err <- max(abs(c(y0err, varerr)))
    expect_true(err < 0.05)
})

test_that("test heuristic expansion on theta=5/4 pi",{
    theta <- 5/4*pi
    Ra <- 10^6
    y0 <- ecd.y0_isomorphic(theta, Ra)
    d <- ecd.polar(Ra, theta)
    y0a <- -2^(2/3)*abs(d@alpha)^(1/3)
    y0err <- abs(y0a/y0-1)
    
    var <- 9/2^(5/3)*abs(d@alpha)^(2/3)
    varerr <- abs(var/ecd.var(d)-1)
    
    expect_true(varerr < 0.05 & y0err < eps)
})

# ---------------------------------------------
