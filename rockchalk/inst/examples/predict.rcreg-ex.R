
library(rockchalk)
dat <- genCorrelatedData(1000, stde=5)

m1 <- lm(y ~ x1 * x2, data=dat)

m1mc <- meanCenter(m1)
summary(m1mc)

m1rc <- residualCenter(m1)
summary(m1rc)


newdf <- apply(dat, 2, summary)
newdf <- as.data.frame(newdf)

predict(m1rc, newdata=newdf)
