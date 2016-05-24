require(glarma)
### Boat Race:
data(OxBoatRace)

y1 <- OxBoatRace$Camwin
n1 <- rep(1, length(OxBoatRace$Year))
Y <- cbind(y1, n1 - y1)
X <- cbind(OxBoatRace$Intercept, OxBoatRace$Diff)
colnames(X) <- c("Intercept", "Weight Diff")

oxcamglm <- glm(Y ~ Diff + I(Diff^2),
                data = OxBoatRace,
                family = binomial(link = "logit"), x = TRUE)

X <- oxcamglm$x

glarmamod <- glarma(Y, X, thetaLags = c(1, 2), type = "Bin", method = "NR",
                    residuals = "Pearson", maxit = 100, grad = 1e-6)

summary(glarmamod)

par(mfrow=c(3,2))
plot(glarmamod)
### model adequate - test using randomized PIT?
rt <- normRandPIT(glarmamod)$rt
par(mfrow = c(2,2))
hist(rt, main = "Histogram of Randomized Residuals",
     xlab = expression(r[t]))
qqnorm(rt, main = "Q-Q Plot of  Randomized Residuals" )
abline(0, 1, lty = 2)
acf(rt, main = "ACF of Randomized Residuals")
pacf(rt, main = "PACF of Randomized Residuals")
