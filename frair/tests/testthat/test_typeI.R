library(frair)
data("gammarus")
context("Testing typeI")

## Using synthetic data
set.seed(1234)
tol <- 0.01 # A multipler on tolerance... Tolerant to 3 s.f.
len <- 60; tim <- 5
datx <- rep(1:len, times=tim)
X1 <- 0.25
X2 <- 0.75
r1 <- rnorm(len*tim, mean = X1, sd = 0.1)
r2 <- rnorm(len*tim, mean = X2, sd = 0.1)
daty1 <- abs(round(r1*datx, 0))
daty2 <- abs(round(r2*datx, 0))
dat <- data.frame(datx,daty1,daty2)

test_that("friar_fit understands a typeI response", {
    m1 <- frair_fit(daty1~datx, data=dat, response='typeI', 
                 start=list(a=0.5), fixed=list(T=1))
    expect_is(m1, 'frfit')
})

test_that("we recover the mean slope of a typeI response", {
    m1 <- frair_fit(daty1~datx, data=dat, response='typeI', 
                    start=list(a=0.5), fixed=list(T=1))
    m2 <- frair_fit(daty2~datx, data=dat, response='typeI', 
                    start=list(a=0.5), fixed=list(T=1))
    a1 <- coef(m1)['a']
    a2 <- coef(m2)['a']
    expect_less_than(abs(a1-X1), X1*tol)
    expect_less_than(abs(a2-X2), X2*tol)
})

test_that("different values for 'T' have the intended outcome", {
    m2 <- frair_fit(daty2~datx, data=dat, response='typeI', 
                    start=list(a=0.5), fixed=list(T=1))
    m3 <- frair_fit(daty2~datx, data=dat, response='typeI', 
                    start=list(a=0.5/2), fixed=list(T = 2))
    m4 <- frair_fit(daty2~datx, data=dat, response='typeI', 
                    start=list(a=0.5/0.5), fixed=list(T = 0.5))
    m5 <- frair_fit(daty2~datx, data=dat, response='typeI', 
                    start=list(a=0.5/48), fixed=list(T = 48))
    a2 <- coef(m2)['a']
    a3 <- coef(m3)['a']*2
    a4 <- coef(m4)['a']*0.5
    a5 <- coef(m5)['a']*48
    expect_less_than(abs(a2-a3), a2*tol)
    expect_less_than(abs(a2-a4), a2*tol)
    expect_less_than(abs(a2-a5), a2*tol)
})

# With real data (albeit an inappropriate model for these data)
## TODO: Gammarus?


