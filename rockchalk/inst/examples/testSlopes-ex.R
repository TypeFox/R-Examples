library(rockchalk)
library(car)
m1 <- lm(statusquo ~ income * age + education + sex + age, data = Chile)
m1ps <- plotSlopes(m1, modx = "income", plotx = "age")
m1psts <- testSlopes(m1ps)
plot(m1psts)


dat2 <- genCorrelatedData(N = 400, rho = .1, means = c(50, -20),
                          stde = 300, beta = c(2, 0, 0.1, -0.4))
m2 <- lm(y ~ x1*x2, data = dat2)
m2ps <- plotSlopes(m2, plotx = "x1", modx = "x2")
m2psts <- testSlopes(m2ps)
plot(m2psts)
m2ps <- plotSlopes(m2, plotx = "x1", modx = "x2", modxVals = "std.dev", n = 5)
m2psts <- testSlopes(m2ps)
plot(m2psts)

## Try again with longer variable names

colnames(dat2) <- c("oxygen","hydrogen","species")
m2a <- lm(species ~ oxygen*hydrogen, data = dat2)
m2aps1 <- plotSlopes(m2a, plotx = "oxygen", modx = "hydrogen")
m2aps1ts <- testSlopes(m2aps1)
plot(m2aps1ts)
m2aps2 <- plotSlopes(m2a, plotx = "oxygen", modx = "hydrogen",
                     modxVals = "std.dev", n = 5)
m2bps2ts <- testSlopes(m2aps2)
plot(m2bps2ts)



dat3 <- genCorrelatedData(N = 400, rho = .1, stde = 300,
                          beta = c(2, 0, 0.3, 0.15),
                          means = c(50,0), sds = c(10, 40))
m3 <- lm(y ~ x1*x2, data = dat3)
m3ps <- plotSlopes(m3, plotx = "x1", modx = "x2")
m3sts <- testSlopes(m3ps)
plot(testSlopes(m3ps))
plot(testSlopes(m3ps), shade = FALSE)

## Finally, if model has no relevant interactions, testSlopes does nothing.
m9 <- lm(statusquo ~ age + income * education + sex + age, data = Chile)
m9ps <- plotSlopes(m9, modx = "education", plotx = "age", plotPoints = FALSE)
m9psts <- testSlopes(m9ps)
