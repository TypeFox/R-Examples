library(rockchalk)

## Replicate some R classics.  The budworm.lg data from predict.glm
## will work properly after re-formatting the information as a data.frame:

## example from Venables and Ripley (2002, pp. 190-2.)
df <- data.frame(ldose = rep(0:5, 2),
                 sex = factor(rep(c("M", "F"), c(6, 6))),
                 SF.numdead = c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16))
df$SF.numalive = 20 - df$SF.numdead

budworm.lg <- glm(cbind(SF.numdead, SF.numalive) ~ sex*ldose,
                  data = df,  family = binomial)

predictOMatic(budworm.lg)

predictOMatic(budworm.lg, n = 7)

predictOMatic(budworm.lg, predVals = c("ldose"), n = 7)

predictOMatic(budworm.lg, predVals = c(ldose = "std.dev.", sex = "table"))



## Now make up a data frame with several numeric and categorical predictors.

set.seed(12345)
N <- 100
x1 <- rpois(N, l = 6)
x2 <- rnorm(N, m = 50, s = 10)
x3 <- rnorm(N)
xcat1 <- gl(2,50, labels = c("M","F"))
xcat2 <- cut(rnorm(N), breaks = c(-Inf, 0, 0.4, 0.9, 1, Inf),
             labels = c("R", "M", "D", "P", "G"))
dat <- data.frame(x1, x2, x3, xcat1, xcat2)
rm(x1, x2, x3, xcat1, xcat2)
dat$xcat1n <- with(dat, contrasts(xcat1)[xcat1, , drop = FALSE])
dat$xcat2n <- with(dat, contrasts(xcat2)[xcat2, ])
STDE <- 15
dat$y <- with(dat,
              0.03 + 0.8*x1 + 0.1*x2 + 0.7*x3 + xcat1n %*% c(2) +
              xcat2n %*% c(0.1,-2,0.3, 0.1) + STDE*rnorm(N))
## Impose some random missings
dat$x1[sample(N, 5)] <- NA
dat$x2[sample(N, 5)] <- NA
dat$x3[sample(N, 5)] <- NA
dat$xcat2[sample(N, 5)] <- NA
dat$xcat1[sample(N, 5)] <- NA
dat$y[sample(N, 5)] <- NA
summarize(dat)


m0 <- lm(y ~ x1 + x2 + xcat1, data = dat)
summary(m0)
## The model.data() function in rockchalk creates as near as possible
## the input data frame.
m0.data <- model.data(m0)
summarize(m0.data)

## no predVals: analyzes each variable separately
(m0.p1 <- predictOMatic(m0))

## requests confidence intervals from the predict function
(m0.p2 <- predictOMatic(m0, interval = "confidence"))

## predVals as vector of variable names: gives "mix and match" predictions
(m0.p3 <- predictOMatic(m0, predVals = c("x1", "x2")))

## predVals as vector of variable names: gives "mix and match" predictions
(m0.p3s <- predictOMatic(m0, predVals = c("x1", "x2"), divider = "std.dev."))

## "seq" is an evenly spaced sequence across the predictor.
(m0.p3q <- predictOMatic(m0, predVals = c("x1", "x2"), divider = "seq"))

(m0.p3i <- predictOMatic(m0, predVals = c("x1", "x2"),
                         interval = "confidence", n = 3))

(m0.p3p <- predictOMatic(m0, predVals = c("x1", "x2"), divider = pretty))

## predVals as vector with named divider algorithms.
(m0.p3 <- predictOMatic(m0, predVals = c(x1 = "seq", x2 = "quantile")))
## predVals as named vector of divider algorithms

## same idea, decided to double-check
(m0.p3 <- predictOMatic(m0, predVals = c(x1 = "quantile", x2 = "std.dev.")))
getFocal(m0.data$x2, xvals =  "std.dev.", n = 5)


## Change from quantile to standard deviation divider
(m0.p5 <- predictOMatic(m0, divider = "std.dev.", n = 5))

## Still can specify particular values if desired
(m0.p6 <- predictOMatic(m0, predVals = list("x1" = c(6,7),
                            "xcat1" = levels(m0.data$xcat1))))

(m0.p7 <- predictOMatic(m0, predVals = c(x1 = "quantile", x2 = "std.dev.")))
getFocal(m0.data$x2, xvals =  "std.dev.", n = 5)

(m0.p8 <- predictOMatic(m0, predVals = list( x1 = quantile(m0.data$x1,
                        na.rm = TRUE, probs = c(0, 0.1, 0.5, 0.8,
                        1.0)), xcat1 = levels(m0.data$xcat1))))

(m0.p9 <- predictOMatic(m0, predVals = list(x1 = "seq", "xcat1" =
                                levels(m0.data$xcat1)), n = 8) )


(m0.p10 <- predictOMatic(m0, predVals = list(x1 = "quantile",
                                 "xcat1" = levels(m0.data$xcat1)), n =  5) )


(m0.p11 <- predictOMatic(m0, predVals = c(x1 = "std.dev."), n = 10))

## Previous same as

(m0.p11 <- predictOMatic(m0, predVals = c(x1 = "default"), divider =
                         "std.dev.", n = 10))

## Previous also same as

(m0.p11 <- predictOMatic(m0, predVals = c("x1"), divider = "std.dev.", n = 10))


(m0.p11 <- predictOMatic(m0,  predVals = list(x1 = c(0, 5, 8), x2 = "default"),
                         divider = "seq"))



m1 <- lm(y ~ log(10+x1) + sin(x2) + x3, data = dat)
m1.data <- model.data(m1)
summarize(m1.data)


(newdata(m1))
(newdata(m1, predVals = list(x1 = c(6, 8, 10))))
(newdata(m1, predVals = list(x1 = c(6, 8, 10), x3 = c(-1,0,1))))
(newdata(m1, predVals = list(x1 = c(6, 8, 10),
                 x2 = quantile(m1.data$x2, na.rm = TRUE), x3 = c(-1,0,1))))

(m1.p1 <- predictOMatic(m1, divider = "std.dev", n = 5))
(m1.p2 <- predictOMatic(m1, divider = "quantile", n = 5))

(m1.p3 <- predictOMatic(m1, predVals = list(x1 = c(6, 8, 10),
                            x2 = median(m1.data$x2, na.rm = TRUE))))

(m1.p4 <- predictOMatic(m1, predVals = list(x1 = c(6, 8, 10),
                                x2 = quantile(m1.data$x2, na.rm = TRUE))))

(m1.p5 <- predictOMatic(m1))
(m1.p6 <- predictOMatic(m1, divider = "std.dev."))
(m1.p7 <- predictOMatic(m1, divider = "std.dev.", n = 3))
(m1.p8 <- predictOMatic(m1, divider = "std.dev.", interval = "confidence"))


m2 <- lm(y ~ x1 + x2 + x3 + xcat1 + xcat2, data = dat)
##  has only columns and rows used in model fit
m2.data <- model.data(m2)
summarize(m2.data)

## Check all the margins
(predictOMatic(m2, interval = "conf"))

## Lets construct predictions the "old fashioned way" for comparison

m2.new1 <- newdata(m2, predVals = list(xcat1 = levels(m2.data$xcat1),
                           xcat2 = levels(m2.data$xcat2)), n = 5)
predict(m2, newdata = m2.new1)


(m2.p1 <- predictOMatic(m2,
                        predVals = list(xcat1 = levels(m2.data$xcat1),
                            xcat2 = levels(m2.data$xcat2)),
                        xcat2 = c("M","D")))
## See? same!

## Pick some particular values for focus
m2.new2 <- newdata(m2, predVals = list(x1 = c(1,2,3), xcat2 = c("M","D")))
## Ask for predictions
predict(m2, newdata = m2.new2)


## Compare: predictOMatic generates a newdata frame and predictions in one step

(m2.p2 <- predictOMatic(m2, predVals = list(x1 = c(1,2,3),
                                xcat2 = c("M","D"))))

(m2.p3 <- predictOMatic(m2, predVals = list(x2 = c(0.25, 1.0),
                                xcat2 = c("M","D"))))

(m2.p4 <- predictOMatic(m2, predVals = list(x2 = plotSeq(m2.data$x2, 10),
                                xcat2 = c("M","D"))))

(m2.p5 <- predictOMatic(m2, predVals = list(x2 = c(0.25, 1.0),
                                xcat2 = c("M","D")), interval = "conf"))

(m2.p6 <- predictOMatic(m2, predVals = list(x2 = c(49, 51),
                                xcat2 = levels(m2.data$xcat2),
                                x1 = plotSeq(dat$x1))))

plot(y ~ x1, data = m2.data)
by(m2.p6, list(m2.p6$xcat2), function(x) {
    lines(x$x1, x$fit, col = x$xcat2, lty = as.numeric(x$xcat2))
})

m2.newdata <- newdata(m2, predVals = list(x2 = c(48, 50, 52),
                              xcat2 = c("M","D")))
predict(m2, newdata = m2.newdata)

(m2.p7 <- predictOMatic(m2, predVals = list(x2 = c(48, 50, 52),
                                xcat2 = c("M","D"))))

(m2.p8 <- predictOMatic(m2,
             predVals = list(x2 = range(m2.data$x2, na.rm = TRUE),
             xcat2 = c("M","D"))))

(m2.p9 <- predictOMatic(m2, predVals = list(x2 = plotSeq(m2.data$x2),
             x1 = quantile(m2.data$x1, pr =c(0.33, 0.66), na.rm = TRUE),
             xcat2 = c("M","D"))))
plot(y ~ x2 , data = m2.data)

by(m2.p9, list(m2.p9$x1, m2.p9$xcat2), function(x) {lines(x$x2, x$fit)})



(predictOMatic(m2, predVals = list(x2 = c(50, 60), xcat2 = c("M","D")),
               interval = "conf"))

## create a dichotomous dependent variable
y2 <- ifelse(rnorm(N) > 0.3, 1, 0)
dat <- cbind(dat, y2)

m3 <- glm(y2 ~ x1 + x2 + x3 + xcat1, data = dat, family = binomial(logit))
summary(m3)
m3.data <- model.data(m3)
summarize(m3.data)

(m3.p1 <- predictOMatic(m3, divider = "std.dev."))

(m3.p2 <- predictOMatic(m3, predVals = list(x2 = c(40, 50, 60),
                             xcat1 = c("M","F")),
                        divider = "std.dev.", interval = "conf"))

## Want a full accounting for each value of x2?
(m3.p3 <- predictOMatic(m3,
                predVals = list(x2 = unique(m3.data$x2),
                    xcat1 = c("M","F")), interval = "conf"))


## Would like to write a more beautiful print method
## for output object, but don't want to obscure structure from user.
## for (i in names(m3.p1)){
##     dns <- cbind(m3.p1[[i]][i], m3.p1[[i]]$fit)
##     colnames(dns) <- c(i, "predicted")
##     print(dns)
## }


