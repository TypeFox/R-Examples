
# convert freeny data set to a zooreg series (without the lag column)
data(freeny)
library(dyn)
freeny.z <- zooreg(cbind(y = as.vector(freeny.y), freeny.x[,-1]),
		start = start(freeny.y), freq = frequency(freeny.y))

freeny.lm <- dyn$lm(y ~ lag(y,-1) + ., freeny.z)
summary(freeny.lm)

y <- freeny.z[,"y"]
y.fit <- fitted(freeny.lm)
plot(merge(y, y.fit), col = 1:2, plot = "single")

library(quantreg)
freeny.rq <- dyn$rq(y ~ lag(y,-1) + ., freeny.z, tau = 1:9/10)
freeny.rq
plot(zoo(t(coef(freeny.rq)), freeny.rq$tau), main = "freeny", xlab = "tau")
head(resid(freeny.rq))

