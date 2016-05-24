### R code from vignette source 'usl.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: options
###################################################
options(prompt="R> ", scipen=4, digits=3, width=66, show.signif.stars=FALSE)
options(SweaveHooks=list(fig = function() par(mar = c(4.1, 4.1, 2.1, 2.1))))
set.seed(42, kind = "default", normal.kind = "default")


###################################################
### code chunk number 2: usl.Rnw:191-194
###################################################
library(usl)
data(raytracer)
raytracer


###################################################
### code chunk number 3: rtplot1
###################################################
plot(throughput ~ processors, data = raytracer)


###################################################
### code chunk number 4: usl.Rnw:216-217
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(throughput ~ processors, data = raytracer)


###################################################
### code chunk number 5: usl.Rnw:231-232
###################################################
usl.model <- usl(throughput ~ processors, data = raytracer)


###################################################
### code chunk number 6: usl.Rnw:237-238
###################################################
summary(usl.model)


###################################################
### code chunk number 7: usl.Rnw:270-271
###################################################
efficiency(usl.model)


###################################################
### code chunk number 8: rtbarplot
###################################################
barplot(efficiency(usl.model), ylab = "efficiency / processor", xlab = "processors")


###################################################
### code chunk number 9: usl.Rnw:284-285
###################################################
getOption("SweaveHooks")[["fig"]]()
barplot(efficiency(usl.model), ylab = "efficiency / processor", xlab = "processors")


###################################################
### code chunk number 10: usl.Rnw:298-299
###################################################
coef(usl.model)


###################################################
### code chunk number 11: usl.Rnw:305-306
###################################################
confint(usl.model, level = 0.95)


###################################################
### code chunk number 12: rtplot2
###################################################
plot(throughput ~ processors, data = raytracer, pch = 16)
plot(usl.model, add = TRUE)


###################################################
### code chunk number 13: usl.Rnw:326-327
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(throughput ~ processors, data = raytracer, pch = 16)
plot(usl.model, add = TRUE)


###################################################
### code chunk number 14: usl.Rnw:339-340
###################################################
predict(usl.model, data.frame(processors = c(96, 128)))


###################################################
### code chunk number 15: usl.Rnw:348-349
###################################################
peak.scalability(usl.model)


###################################################
### code chunk number 16: usl.Rnw:368-371
###################################################
library(usl)
data(specsdm91)
specsdm91


###################################################
### code chunk number 17: usl.Rnw:384-385
###################################################
usl.model <- usl(throughput ~ load, specsdm91, method = "nlxb")


###################################################
### code chunk number 18: usl.Rnw:415-416
###################################################
summary(usl.model)


###################################################
### code chunk number 19: usl.Rnw:434-436
###################################################
peak.scalability(usl.model)
peak.scalability(usl.model, kappa = 0.00005)


###################################################
### code chunk number 20: spplot1
###################################################
plot(specsdm91, pch = 16, ylim = c(0,2500))
plot(usl.model, add = TRUE)
cache.scale <- scalability(usl.model, kappa = 0.00005)
curve(cache.scale, lty = 2, add = TRUE)


###################################################
### code chunk number 21: usl.Rnw:474-475
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(specsdm91, pch = 16, ylim = c(0,2500))
plot(usl.model, add = TRUE)
cache.scale <- scalability(usl.model, kappa = 0.00005)
curve(cache.scale, lty = 2, add = TRUE)


###################################################
### code chunk number 22: usl.Rnw:489-492
###################################################
scalability(usl.model)(peak.scalability(usl.model))
scf <- scalability(usl.model, kappa = 0.00005)
scf(peak.scalability(usl.model, kappa = 0.00005))


###################################################
### code chunk number 23: usl.Rnw:504-505
###################################################
load <- with(specsdm91, expand.grid(load = seq(min(load), max(load))))


###################################################
### code chunk number 24: usl.Rnw:513-514
###################################################
fit <- predict(usl.model, newdata = load, interval = "confidence", level = 0.95)


###################################################
### code chunk number 25: usl.Rnw:522-525
###################################################
usl.polygon <- matrix(c(load[, 1], rev(load[, 1]),
                        fit[, 'lwr'], rev(fit[, 'upr'])),
                 nrow = 2 * nrow(load))


###################################################
### code chunk number 26: ciplot1
###################################################
# Create empty plot (define canvas size, axis, ...)
plot(specsdm91,
     ylim = c(0, 2000), 
     xlab = names(specsdm91)[1],
     ylab = names(specsdm91)[2], 
     type = "n")

# Plot gray polygon indicating the confidence interval
polygon(usl.polygon, border = NA, col = "gray")

# Plot the measured throughput
points(specsdm91, pch = 16)

# Plot the fit
lines(load[, 1], fit[, 'fit'])


###################################################
### code chunk number 27: usl.Rnw:555-556
###################################################
getOption("SweaveHooks")[["fig"]]()
# Create empty plot (define canvas size, axis, ...)
plot(specsdm91,
     ylim = c(0, 2000), 
     xlab = names(specsdm91)[1],
     ylab = names(specsdm91)[2], 
     type = "n")

# Plot gray polygon indicating the confidence interval
polygon(usl.polygon, border = NA, col = "gray")

# Plot the measured throughput
points(specsdm91, pch = 16)

# Plot the fit
lines(load[, 1], fit[, 'fit'])


###################################################
### code chunk number 28: usl.Rnw:599-602
###################################################
load <- data.frame(load = c(10, 20, 100, 200))
ovhd <- overhead(usl.model, newdata = load)
ovhd


###################################################
### code chunk number 29: ovplot1
###################################################
barplot(height = t(ovhd), names.arg = load[, 1],
        xlab = names(load), legend.text = TRUE)


###################################################
### code chunk number 30: usl.Rnw:631-632
###################################################
getOption("SweaveHooks")[["fig"]]()
barplot(height = t(ovhd), names.arg = load[, 1],
        xlab = names(load), legend.text = TRUE)


