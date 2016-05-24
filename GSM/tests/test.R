library(GSM)

# simulate some data
set.seed(2040)
y <- rgsm(500, c(.1, .3, .4, .2), 1)

# fit the model
burnin <- 500
mcmcsim <- 1000
J <- 300
gsm.out <- estim.gsm(y, J, 300, burnin + mcmcsim, 6500, 340, 1/J)

# summary of fit
summary(gsm.out, TRUE, start = (burnin + 1))

# plot of fit
par(mfrow = c(3, 2))
plot(gsm.out)
plot(gsm.out, ndens = 0, nbin = 20, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, histogram = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, histogram = TRUE, bands = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 5, nbin = 20, histogram = TRUE, bands = TRUE, start = (burnin + 1))
plot(gsm.out, ndens = 0, nbin = 20, bands = TRUE, start = (burnin + 1))

# estimation of the tail probabilities
par(mfrow = c(1, 1))
thresh <- c(0.1, 0.5, 0.75, 1, 2)
tail.prob.est <- tail.prob.true <- rep(NA, length(thresh))
for (i in 1:length(thresh)){
   tail.prob.est[i] <- mean(predict(gsm.out, thresh[i]))
   tail.prob.true[i] <- sum(y > thresh[i])/length(y)
}
qqplot(tail.prob.true, tail.prob.est, main = "Q-Q plot of true vs. estimated tail probability")
abline(0, 1, lty = 2)
