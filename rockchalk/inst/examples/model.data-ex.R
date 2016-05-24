library(rockchalk)






## first, check if model.data works when there is no data argument
## This used to fail, now OK

x1 <- rnorm(100, m = 100, s = 10)
x2 <- rnorm(100, m = 50, s =20)
y <- rnorm(100, m = 40, s = 3)

m0 <- lm(y ~ log(10+x1) + x2)
m0.data <- model.data(m0)
head(m0.data)


m1 <- lm(log(43 + y) ~ log(10+x1) + x2)
m1.data <- model.data(m1)
head(m1.data)

d <- 3

m2 <- lm(log(d + y) ~ log(10+x1) + x2)
m2.data <- model.data(m2)
head(m2.data)

m3 <- lm(log(y + d) ~ log(10+x1) + x2)
m3.data <- model.data(m3)
head(m3.data)



## check numeric and categorical predictors

x1 <- rpois(100, l=6)
x2 <- rnorm(100, m=50, s=10)
x3 <- rnorm(100)
xcat1 <- gl(2,50, labels=c("M","F"))
xcat2 <- cut(rnorm(100), breaks=c(-Inf, 0, 0.4, 0.9, 1, Inf),
             labels=c("R", "M", "D", "P", "G"))
dat <- data.frame(x1, x2, x3, xcat1, xcat2)
rm(x1, x2, x3, xcat1, xcat2)
dat$xcat1n <- with(dat, contrasts(xcat1)[xcat1, ,drop=FALSE])
dat$xcat2n <- with(dat, contrasts(xcat2)[xcat2, ])

STDE <- 20
dat$y <- with(dat,
              0.03 + 0.8*x1 + 0.1*x2 + 0.7*x3 +
              xcat1n %*% c(2) + xcat2n %*% c(0.1,-2,0.3, 0.1) +
              STDE*rnorm(100))



m1 <- lm(y ~ poly(x1, 2), data=dat)
m1.data <- model.data(m1)
head(m1.data)
attr(m1.data, "varNamesRHS")

## Check to make sure d is not mistaken for a data column
d <- 2
m2 <- lm(y ~ poly(x1, d), data=dat)
m2.data <- model.data(m2)
head(m2.data)
attr(m2.data, "varNamesRHS")


## Check to see how the 10 in log is handled
m3 <- lm(y ~ log(10 + x1) + poly(x1, d) + sin(x2), data=dat)
m3.data <- model.data(m3)
head(m3.data)
attr(m3.data, "varNamesRHS")


m4 <- lm(log(50+y) ~ log(d+10+x1) + poly(x1, 2), data=dat)
m4.data <- model.data(m4)
head(m4.data)
attr(m4.data, "varNamesRHS")


m5 <- lm(y ~ x1*x1, data=dat)
m5.data <- model.data(m5)
head(m5.data)
attr(m5.data, "varNamesRHS")


m6 <- lm(y ~ x1 + I(x1^2), data=dat)
m6.data <- model.data(m6)
head(m6.data)
attr(m6.data, "varNamesRHS")


## Put in some missings.
## poly doesn't work if there are missings, but
## can test with log
dat$x1[sample(100, 5)] <- NA
dat$y[sample(100, 5)] <- NA
dat$x2[sample(100, 5)] <- NA
dat$x3[sample(100,10)] <- NA

m1 <- lm(y ~ log(10 + x1), data=dat)
m1.data <- model.data(m1)
head(m1.data)
summarize(m1.data)
attr(m1.data, "varNamesRHS")



m2 <- lm(y ~ log(x1 + 10), data=dat)
m2.data <- model.data(m2)
head(m2.data)
summarize(m1.data)
attr(m1.data, "varNamesRHS")

d <- 2
m3 <- lm(log(50+y) ~ log(d+10+x1) + x2 + sin(x3), data=dat)
m3.data <- model.data(m3)
head(m3.data)
summarize(m3.data)
attr(m3.data, "varNamesRHS")


m4 <- lm(y ~ I(x1) + I(x1^2) + log(x2), data=dat)
m4.data <- model.data(m4)
summarize(m4.data)
attr(m4.data, "varNamesRHS")


m5 <- lm(y ~ x1 + I(x1^2) + cos(x2), data=dat)
m5.data <- model.data(m5)
head(m5.data)
summarize(m5.data)
attr(m5.data, "varNamesRHS")



## Now try with some variables in the dataframe, some not

x10 <- rnorm(100)
x11 <- rnorm(100)



m6 <- lm(y ~ x1 + I(x1^2) + cos(x2) + log(10 + x10) + sin(x11) +
         x10*x11, data = dat)
m6.data <- model.data(m6)
head(m6.data)
dim(m6.data)
summarize(m5.data)
attr(m6.data, "varNamesRHS")


