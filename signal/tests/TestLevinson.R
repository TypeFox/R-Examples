library("signal")
load("savedTestLev.Rdata")

set.seed(234)
y <- arima.sim(n = 1000, list(ar = c(1, -0.9, 0.3)))
z <- arima.sim(n = 1000, list(ar = c(1, -0.4, 0.1)))
x <- cbind(as.numeric(y), as.numeric(z))
TestLev1 <- levinson(x, 2)
TestLev2 <- levinson(y, 2)
TestLev3 <- levinson(z, 2)

all.equal(savedTestLev1, TestLev1)
all.equal(savedTestLev2, TestLev2)
all.equal(savedTestLev3, TestLev3)

