### R code from vignette source 'UserGuide.Rnw'

###################################################
### code chunk number 1: UserGuide.Rnw:42-44
###################################################
library(schwartz97)
options(width = 60)  


###################################################
### code chunk number 2: UserGuide.Rnw:239-252
###################################################
s0 <- 100
delta0 <- 0
mu <- 0.1
sigmaS <- 0.2
kappa <- 1
alpha <- 0.1
sigmaE <- 0.3
rho <- 0.4

obj <- schwartz2f(s0 = s0, delta0 = delta0, alpha = alpha,
                       mu = mu, sigmaS = sigmaS, sigmaE = sigmaE,
                       rho = rho, kappa = kappa)
obj


###################################################
### code chunk number 3: UserGuide.Rnw:277-282
###################################################
time <- 5
sample.t <- rstate(n = 2000, time, obj)
pstate(c(0, -Inf), c(150, 0), time, obj)
mean(obj, time = c(1, 10))
plot(obj, n = 30, time = 5, dt = 1 / 52)


###################################################
### code chunk number 4: UserGuide.Rnw:286-287
###################################################
plot(obj, n = 50, time = 5, dt = 1 / 52)


###################################################
### code chunk number 5: UserGuide.Rnw:323-339
###################################################
s0 <- 80
delta0 <- 0.05
mu <- 0.1
sigmaS <- 0.3
kappa <- 1.5
alpha <- 0.05
sigmaE <- 0.4
rho <- 0.6
lambda <- 0.04
r <- 0.03
set.seed(1)
obj <- schwartz2f(s0, delta0, mu, sigmaS, kappa, alpha, sigmaE, rho)
state.traj <- simstate(n = 52 * time, time, obj)
pricefutures(seq(0, 2, by = 0.4), obj, lambda = lambda, r = r)
priceoption(type = "call", time = 1, Time = 2, K = 85,
            obj, r = r, lambda = lambda)


###################################################
### code chunk number 6: UserGuide.Rnw:343-357
###################################################
ttm <- seq(0, 2, by = 0.4)
forward.curves <- sapply(ttm, pricefutures, s0 = state.traj[,1], delta0 = state.traj[,2],
                         sigmaS = sigmaS, kappa = kappa, alpha = alpha, sigmaE = sigmaE,
                         lambda = lambda, r = r)
ttm.mat <- matrix(ttm, ncol = 6, nrow = 52 * time, byrow = TRUE)
dates <- Sys.Date() + seq(0, length = time * 52, by = 7)
rownames(ttm.mat) <- as.character(dates)
par(oma = c(6, 0, 5, 0), mfrow = c(2, 1), mar = c(0, 5.1, 0, 2.1))
## Prepare a list which is accepted by 'futuresplot'
futuresplot(list(price = forward.curves, ttm = ttm.mat * 260), xaxt = "n")
plot(dates, state.traj[,2], main = "", type = "l", ylab = "Convenience yield",
     xlim = c(dates[1], rev(dates)[1] + max(ttm) * 260))
abline(h = 0)



###################################################
### code chunk number 7: UserGuide.Rnw:389-413
###################################################
s0 <- 1
delta0 <- 0.0
sigmaS <- 0.3
kappa <- 1
sigmaE <- 0.4
rho <- 0.5
r <- 0.03
ttm <- 0:4
## Pure contango
pricefutures(ttm, s0 = s0, delta0 = 0, sigmaS = sigmaS,
             kappa = kappa, sigmaE = sigmaE, rho = rho,
             r = r, alphaT = 0)
## Backwardation and then contango
pricefutures(ttm, s0 = s0, delta0 = 2 * r, sigmaS = sigmaS,
             kappa = kappa, sigmaE = sigmaE, rho = rho,
             r = r, alphaT = 0)
## Pure backwardation
pricefutures(ttm, s0 = s0, delta0 = r, sigmaS = sigmaS,
             kappa = kappa, sigmaE = sigmaE, rho = rho,
             r = r, alphaT = 2 * r)
## Contango and then backwardation
pricefutures(ttm, s0 = s0, delta0 = -r, sigmaS = sigmaS,
             kappa = kappa, sigmaE = sigmaE, rho = rho,
             r = r, alphaT = 2 * r)


###################################################
### code chunk number 8: UserGuide.Rnw:417-449
###################################################
s0 <- 1
delta0 <- 0.0
sigmaS <- 0.3
kappa <- 1
sigmaE <- 0.4
rho <- 0.5
r <- 0.03
ttm <- seq(0, 4, by = 0.05)
## Pure contango
plot(ttm, ttm, type = "n", ylim = c(0.9, 1.18), xlab = "Time to maturity [y]",
     ylab = "Futures price", main = "Term structure shapes")
lines(ttm, pricefutures(ttm, s0 = s0, delta0 = 0, sigmaS = sigmaS,
                        kappa = kappa, sigmaE = sigmaE, rho = rho,
                        lambda = 0, r = r, alpha = 0), col = "blue", lty = 2)
## Backwardation and then contango
lines(ttm,pricefutures(ttm, s0 = s0, delta0 = 2 * r, sigmaS = sigmaS,
                       kappa = kappa, sigmaE = sigmaE, rho = rho,
                       lambda = 0, r = r, alpha = 0), col = "blue")
## Pure backwardation
lines(ttm,pricefutures(ttm, s0 = s0, delta0 = r, sigmaS = sigmaS,
                       kappa = kappa, sigmaE = sigmaE, rho = rho,
                       lambda = 0, r = r, alpha = 2 * r), lty = 2)
## Contango and then backwardation
lines(ttm,pricefutures(ttm, s0 = s0, delta0 = -r, sigmaS = sigmaS,
                       kappa = kappa, sigmaE = sigmaE, rho = rho,
                       lambda = 0, r = r, alpha = 2 * r))
legend("topleft", c("delta0 = 0%, alpha = 0%",
                    "delta0 = 6%, alpha = 0%",
                    "delta0 = 3%, alpha = 6%",
                    "delta0 = -3%, alpha = 6%"),
       title = "Common parameters: s0 = 1, sigmaS = 30%, kappa = 1, sigmaE = 40%, rho = 50%, lambda = 0, r = 3%",
       lty = c(2, 1, 2, 1), col = c(rep("blue", 2), rep("black", 2)), lwd = 2, bty = "n")


###################################################
### code chunk number 9: UserGuide.Rnw:476-477
###################################################
  args(fit.schwartz2f)


###################################################
### code chunk number 10: UserGuide.Rnw:563-568
###################################################
data(futures)
wheat.fit <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
                          deltat = 1 / 260, control = list(maxit = 300), silent = TRUE)
wheat.fit
plot(wheat.fit, type = "trace.pars")


###################################################
### code chunk number 11: UserGuide.Rnw:572-573
###################################################
plot(wheat.fit, type = "trace.pars")


###################################################
### code chunk number 12: UserGuide.Rnw:601-612
###################################################
vol.std <- colSums(futures$wheat$vol, na.rm = TRUE) / sum(futures$wheat$vol, na.rm = TRUE)
wheat.fit.constr <- fit.schwartz2f(futures$wheat$price, futures$wheat$ttm / 260,
                                   kappa = 1,
                                   opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE, 
                                     sigmaS = TRUE, kappa = FALSE, alpha = TRUE, 
                                     sigmaE = TRUE, rho = TRUE, lambda = FALSE),
                                   meas.sd = 1 / vol.std  / sum(1 / vol.std) * length(vol.std) * 0.01,
                                   deltat = 1 / 260, control = list(maxit = 300), silent = TRUE)

wheat.fit.constr
plot(wheat.fit.constr, type = "trace.pars")


###################################################
### code chunk number 13: UserGuide.Rnw:615-616
###################################################
plot(wheat.fit.constr, type = "trace.pars")


###################################################
### code chunk number 14: UserGuide.Rnw:634-635
###################################################
plot(wheat.fit.constr, type = "sim", time = 5, n = 30)


###################################################
### code chunk number 15: UserGuide.Rnw:643-652
###################################################
wheat.2007 <- lapply(futures$wheat,
                     function(x)x[as.Date(rownames(x)) > "2007-01-01" & as.Date(rownames(x)) < "2008-07-01",])
par(mfrow = c(1, 2))
futuresplot(wheat.2007, type = "forward.curve")
plot(wheat.fit.constr, type = "forward.curve", data = wheat.2007$price,
     ttm = wheat.2007$ttm / 260)
##xx <- filter.schwartz2f(data = wheat.2007$price,
##                        ttm = wheat.2007$ttm / 260, wheat.fit.constr)
##plot(ts(xx$state, frequency = 260, start = 2007))


###################################################
### code chunk number 16: UserGuide.Rnw:656-660
###################################################
par(mfrow = c(1, 2))
futuresplot(wheat.2007, type = "forward.curve")
plot(wheat.fit.constr, type = "forward.curve", data = wheat.2007$price,
     ttm = wheat.2007$ttm / 260)


###################################################
### code chunk number 17: UserGuide.Rnw:695-701
###################################################
model.resid <- resid(wheat.fit.constr, data = futures$wheat$price, ttm = futures$wheat$ttm / 260,
                     type = "filter.std")
acf(model.resid, na.action = na.pass)

par(mfrow = c(3, 2))
invisible(apply(model.resid, 2, function(x)plot(density(na.omit(x)))))


###################################################
### code chunk number 18: UserGuide.Rnw:705-706
###################################################
acf(model.resid, na.action = na.pass)


###################################################
### code chunk number 19: UserGuide.Rnw:722-724
###################################################
par(mfrow = c(3, 2))
dummy <- apply(model.resid, 2, function(x)plot(density(na.omit(x))))


###################################################
### code chunk number 20: UserGuide.Rnw:743-769
###################################################
state <- filter.schwartz2f(data = futures$wheat$price, ttm = futures$wheat$ttm / 260, wheat.fit.constr)$state
coefs <- coef(wheat.fit.constr)
n <- nrow(futures$wheat$price)
q.fut <- sapply(futures$wheat$ttm[n,] / 260, function(ttm, ...)
                qfutures(ttm = ttm, ...),
                p = c(0.05, 0.95), time = 5 / 260, s0 = state[n,1], delta0 = state[n,2],
                mu = coefs$mu, sigmaS = coefs$sigmaS, kappa = coefs$kappa, alpha = coefs$alpha,
                sigmaE = coefs$sigmaE, rho = coefs$rho, r = coefs$r, lambda = coefs$lambda)

plot(futures$wheat$ttm[n,], futures$wheat$price[n,], ylim = c(650, 850), type = "b",
     xlab = "Time to maturity [d]", ylab = "Price")
points(futures$wheat$ttm[n,], q.fut[1,], col = "blue", type = "b")
points(futures$wheat$ttm[n,], q.fut[2,], col = "blue", type = "b")
legend("topleft", c("Current observed futures price", "One week ahead 90% confidence interval"),
       fill = c("black", "blue"))
## q.fut <- sapply(1:nrow(state), function(i, s0, delta0, ttm, ...)
##                 qfutures(s0 = s0[i], delta0 = delta0[i], ttm = ttm[i], ...),
##                 p = 0.05, time = 1 / 12, ttm = futures$wheat$ttm / 260, s0 = state[,1],
##                 delta0 = state[,2], mu = coefs$mu, sigmaS = coefs$sigmaS, kappa = coefs$kappa, alpha = coefs$alpha,
##                 sigmaE = coefs$sigmaE, rho = coefs$rho, r = coefs$r, lambda = coefs$lambda)

## par(mfrow = c(1, 2))
## plot(ts(cbind(state[,1], futures$wheat$price[,1], q.fut), freq = 260, start = 1995),
##      plot.type = "single", col = c("black", "blue", "red"))
## plot(ts(state[,2], freq = 260, start = 1995))
## abline(h = coefs$alphaT)


###################################################
### code chunk number 21: UserGuide.Rnw:773-779
###################################################
plot(futures$wheat$ttm[n,], futures$wheat$price[n,], ylim = c(650, 850), type = "b",
     xlab = "Time to maturity [d]", ylab = "Price")
points(futures$wheat$ttm[n,], q.fut[1,], col = "blue", type = "b")
points(futures$wheat$ttm[n,], q.fut[2,], col = "blue", type = "b")
legend("topleft", c("Current observed futures price", "One week ahead 90% confidence interval"),
       col = c("black", "blue"), lty = 1)


###################################################
### code chunk number 22: UserGuide.Rnw:805-821
###################################################
futures.w <- rapply(futures, function(x)x[format(as.Date(rownames(x)), "%w") == 3,],
                    classes = "matrix", how = "list")
soybean.meal.fit <- fit.schwartz2f(data = futures.w$soybean.meal$price,
                                   ttm = futures.w$soybean.meal$ttm / 260,
                                   kappa = 1, mu = 0,
                                   opt.pars = c(s0 = FALSE, delta0 = FALSE, mu = TRUE, 
                                     sigmaS = TRUE, kappa = FALSE, alpha = TRUE, 
                                     sigmaE = TRUE, rho = TRUE, lambda = FALSE),
                                   opt.meas.sd = "all", deltat = 1 / 52,
                                   control = list(maxit = 1000), silent = TRUE)
soybean.meal.fit

par(mfrow = c(1, 2))
futuresplot(futures.w$soybean.meal, type = "forward.curve")
plot(soybean.meal.fit, type = "forward.curve", data = futures.w$soybean.meal$price,
     ttm = futures.w$soybean.meal$ttm / 260)


###################################################
### code chunk number 23: UserGuide.Rnw:825-829
###################################################
par(mfrow = c(1, 2))
futuresplot(futures.w$soybean.meal, type = "forward.curve")
plot(soybean.meal.fit, type = "forward.curve", data = futures.w$soybean.meal$price,
     ttm = futures.w$soybean.meal$ttm / 260)


