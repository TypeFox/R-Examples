library(rockchalk)


set.seed(12345)
x1 <- rnorm(100)
x2 <- rnorm(100)
x3 <- rnorm(100)
x4 <- rnorm(100)
y <- rnorm(100)
y2 <- 0.03 + 0.1*x1 + 0.1*x2 + 0.25*x1*x2 + 0.4*x3 -0.1*x4 + 1*rnorm(100)
dat <- data.frame(x1,x2,x3,x4,y, y2)
rm(x1, x2, x3, x4, y, y2)

## linear ordinary regression
m1 <- lm(y ~ x1 + x2 +x3 + x4, data = dat)

plotPlane(m1, plotx1 = "x3", plotx2 = "x4")

plotPlane(m1, plotx1 = "x3", plotx2 = "x4", drawArrows = TRUE)

plotPlane(m1, plotx1 = "x1", plotx2 = "x4", drawArrows = TRUE)


plotPlane(m1, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE, npp = 10)
plotPlane(m1, plotx1 = "x3", plotx2 = "x2", drawArrows = TRUE, npp = 40)

plotPlane(m1, plotx1 = "x3", plotx2 = "x2", drawArrows = FALSE,
          npp = 5, ticktype = "detailed")


## regression with interaction
m2 <- lm(y ~ x1  * x2 +x3 + x4, data = dat)

plotPlane(m2, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE)


plotPlane(m2, plotx1 = "x1", plotx2 = "x4", drawArrows = TRUE)
plotPlane(m2, plotx1 = "x1", plotx2 = "x3", drawArrows = TRUE)

plotPlane(m2, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE,
          phi = 10, theta = 30)



## regression with quadratic;
## Required some fancy footwork in plotPlane, so be happy
dat$y3 <- 0 + 1 * dat$x1 + 2 * dat$x1^2 + 1 * dat$x2 +
    0.4*dat$x3 + 8 * rnorm(100)
m3 <- lm(y3 ~ poly(x1,2) + x2 +x3 + x4, data = dat)
summary(m3)

plotPlane(m3, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE,
          x1lab = "my great predictor", x2lab = "a so-so predictor",
          ylab = "Most awesomest DV ever")

plotPlane(m3, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE,
          x1lab = "my great predictor", x2lab = "a so-so predictor",
          ylab = "Most awesomest DV ever", phi = -20)

plotPlane(m3, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE,
          phi = 10, theta = 30)

plotPlane(m3, plotx1 = "x1", plotx2 = "x4", drawArrows = TRUE,
          ticktype = "detailed")
plotPlane(m3, plotx1 = "x1", plotx2 = "x3", drawArrows = TRUE)

plotPlane(m3, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE,
          phi = 10, theta = 30)

m4 <- lm(y ~ sin(x1) + x2*x3 +x3 + x4, data = dat)
summary(m4)


plotPlane(m4, plotx1 = "x1", plotx2 = "x2", drawArrows = TRUE)
plotPlane(m4, plotx1 = "x1", plotx2 = "x3", drawArrows = TRUE)



eta3 <- 1.1 + .9*dat$x1 - .6*dat$x2 + .5*dat$x3
dat$y4 <- rbinom(100, size = 1, prob = exp( eta3)/(1+exp(eta3)))
gm1 <- glm(y4 ~ x1 + x2 + x3, data = dat, family = binomial(logit))
summary(gm1)
plotPlane(gm1, plotx1 = "x1", plotx2 = "x2")
plotPlane(gm1, plotx1 = "x1", plotx2 = "x2", phi = -10)

plotPlane(gm1, plotx1 = "x1", plotx2 = "x2", ticktype = "detailed")
plotPlane(gm1, plotx1 = "x1", plotx2 = "x2", ticktype = "detailed",
          npp = 30, theta = 30)
plotPlane(gm1, plotx1 = "x1", plotx2 = "x3", ticktype = "detailed",
          npp = 70, theta = 60)

plotPlane(gm1, plotx1 = "x1", plotx2 = "x2", ticktype = c("detailed"),
          npp = 50, theta = 40)

dat$x2 <- 5 * dat$x2
dat$x4 <- 10 * dat$x4
eta4 <- 0.1 + .15*dat$x1 - 0.1*dat$x2 + .25*dat$x3 + 0.1*dat$x4
dat$y4 <- rbinom(100, size = 1, prob = exp( eta4)/(1+exp(eta4)))
gm2 <- glm(y4 ~ x1 + x2 + x3 + x4, data = dat, family = binomial(logit))
summary(gm2)
plotPlane(gm2, plotx1 = "x1", plotx2 = "x2")
plotPlane(gm2, plotx1 = "x2", plotx2 = "x1")
plotPlane(gm2, plotx1 = "x1", plotx2 = "x2", phi = -10)
plotPlane(gm2, plotx1 = "x1", plotx2 = "x2", phi = 5, theta = 70, npp = 40)

plotPlane(gm2, plotx1 = "x1", plotx2 = "x2", ticktype = "detailed")
plotPlane(gm2, plotx1 = "x1", plotx2 = "x2", ticktype = "detailed",
          npp = 30, theta = -30)
plotPlane(gm2, plotx1 = "x1", plotx2 = "x3", ticktype = "detailed",
          npp = 70, theta = 60)

plotPlane(gm2, plotx1 = "x4", plotx2 = "x3", ticktype = "detailed",
          npp = 50, theta = 10)

plotPlane(gm2, plotx1 = "x1", plotx2 = "x2", ticktype = c("detailed"))
