suppressMessages(library(cubfits, quietly = TRUE))
library(EMCluster)
set.seed(1234)

### Get individual of phi.Obs.
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs <- GM / sum(GM) * 15000

### Run ASL optimization.
X <- log(as.matrix(phi.Obs))
ret <- asl.optim(X)
logL.asl <- sum(dasl(X, theta = ret[1], mu = ret[2], sigma = ret[4],
                     log = TRUE))
AIC.asl <- -2 * logL.asl + 2 * 3
BIC.asl <- -2 * logL.asl + log(length(X)) * 3
cat("ASL, logL = ", logL.asl,
    ", AIC = ", AIC.asl, ", BIC = ", BIC.asl, "\n", sep = "")

### Find lognormal/normal.
mu.norm <- mean(X)
sd.norm <- sd(X)
logL.norm <- sum(dnorm(X, mu.norm, sd = sd.norm, log = TRUE))
AIC.norm <- -2 * logL.norm + 2 * 2
BIC.norm <- -2 * logL.norm + log(length(X)) * 2
cat("Normal, logL = ", logL.norm,
    ", AIC = ", AIC.norm, ", BIC = ", BIC.norm, "\n", sep = "")

### GM averaged
x <- seq(min(X), max(X), length = 100)
d.asl <- lapply(x, dasl, theta = ret[1], mu = ret[2], sigma = ret[4])
d.normal <- dnorm(x, mu.norm, sd = sd.norm)

ylim <- c(0, 0.5)
hist(X, nclass = 50, freq = FALSE, ylim = ylim,
     xlab = "Expression (log)", ylab = "Density",
     main = paste("Yassour (N = ", length(X), ", averaged)",
                  sep = ""))
lines(x = x, y = d.normal, col = 1, lty = 1, lwd = 1.5)
lines(x = x, y = d.asl, col = 2, lty = 2, lwd = 1.5)
legend(min(x) + 0.5, ylim[2] * 0.9,
       c("Normal", "ASL"),
       col = 1:2, lty = 1:2, lwd = 1.5)
