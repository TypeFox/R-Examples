library(tsallisqexp)

# ?EPD

#####
# (1) density function
x <- seq(-5, 5, length=101)

cbind(x, y <- dtsal(x, 1/2, 1/4), dtsal(x, 1/2, 1/4, log=TRUE))
# plot(x, y, type="l")
cbind(x, y <- dtsal.tail(x, 1/2, 1/4, xmin=3), dtsal.tail(x, 1/2, 1/4, log=TRUE, xmin=3))

#####
# (2) distribution function

ptsal(x, 1/2, 1/4)
ptsal(x, 1/2, 1/4, lower=FALSE)
ptsal(x, 1/2, 1/4, log=TRUE)

ptsal(x, q=1/2, kappa=4)

ptsal.tail(x, 1/2, 1/4, xmin=3)
ptsal.tail(x, 1/2, 1/4, xmin=3, log=TRUE)
ptsal.tail(x, 1/2, 1/4, xmin=3, lower=FALSE)
ptsal.tail(x, 1/2, 1/4, xmin=3, lower=FALSE, log=TRUE)



#####
# (3) quantile function

qtsal(0:10/10, 3, 2)
qtsal(log(0:10/10), 3, 2, log=TRUE)

qtsal.tail(0:10/10, 3, 2, xmin=3)
qtsal.tail(log(0:10/10), 3, 2, xmin=3, log=TRUE)


#####
# (4) random generation function

rtsal(10, 3, 2)
rtsal.tail(10, 3, 2, xmin=3)

#####
# (5) fit function

set.seed(1234)
x <- rtsal(10, 3, 2)

tsal.fit(x, method="mle.equation")
tsal.fit(x, method="mle.direct")
tsal.fit(x, method="leastsquares")



#####
# (6) boot functions

# ?tsal.boot

tsal.bootstrap.errors(dist=NULL, reps=100, confidence=0.95, n=10)

tsal.bootstrap.errors(dist=tsal.fit(x, method="mle.equation"), reps=100)

tsal.total.magnitude(dist=NULL, n=10)

tsal.total.magnitude(dist=tsal.fit(x, method="mle.equation"))


#####
# (7) test functions

# ?tsal.test

test.tsal.quantile.transform(from=0, to=1e6, shape=1, scale=1,
    n=1e5, lwd=0.01, xmin=0)

test.tsal.LR.distribution(n=10, reps=100, shape=2, scale=3/2,
    xmin=0,method="mle.equation")
test.tsal.LR.distribution(n=1000, reps=100, shape=2, scale=3/2,
    xmin=0,method="mle.equation")
