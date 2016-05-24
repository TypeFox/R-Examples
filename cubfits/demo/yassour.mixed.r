suppressMessages(library(cubfits, quietly = TRUE))
set.seed(1234)

### Get individual of phi.Obs.
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs.all <- yassour[, -1] / sum(GM) * 15000
phi.Obs.all[phi.Obs.all == 0] <- NA

### Get geometric mean of phi.Obs and scaling similar to Wallace (2013).
GM <- apply(yassour[, -1], 1, function(x) exp(mean(log(x[x != 0]))))
phi.Obs <- GM / sum(GM) * 15000

### Run optimization.
X <- log(as.matrix(phi.Obs.all))
param.init <- list(K = 2, prop = c(0.95, 0.05), mu = c(-0.59, 3.11),
                   sigma2 = c(1.40, 0.59), sigma2.e = 0.03)
ret <- mixnormerr.optim(X, K = 2, param = param.init)
print(ret)

### Compute density from 2 mixture normal.
x <- seq(min(X, na.rm = TRUE), max(X, na.rm = TRUE), length = 100)
d.mn <- dmixnormerr(x, ret$param)

### Compute density from posterior of mean expression.
d.post.fits <- dnorm(x, -0.4848434, 0.8731349)
d.post.appr <- dnorm(x, -0.6267638, 1.010109)

### Ploting
par(mfrow = c(2, 1))
ylim <- range(c(0, d.mn, d.post.fits, d.post.appr))

### All data
hist(X, nclass = 100, freq = FALSE, ylim = ylim,
     xlab = "Expression (log)", ylab = "Density",
     main = paste("Yassour (N = ", nrow(X), ", R = ", ncol(X), ")",
                  sep = ""))
lines(x = x, y = d.mn, col = 1, lwd = 1.5)
lines(x = x, y = d.post.fits, col = 2, lty = 2, lwd = 1.5)
lines(x = x, y = d.post.appr, col = 4, lty = 3, lwd = 1.5)
legend(min(x) + 0.5, ylim[2] * 0.9, c("mixture", "fits", "appr"),
       col = c(1, 2, 4), lty = c(1, 2, 3), lwd = 1.5)

### GM averaged
X <- log(phi.Obs)
hist(X, nclass = 100, freq = FALSE, ylim = ylim,
     xlab = "Expression (log)", ylab = "Density",
     main = paste("Yassour (N = ", length(X), ", averaged)",
                  sep = ""))
lines(x = x, y = d.mn, col = 1, lwd = 1.5)
lines(x = x, y = d.post.fits, col = 2, lty = 2, lwd = 1.5)
lines(x = x, y = d.post.appr, col = 4, lty = 3, lwd = 1.5)
legend(min(x) + 0.5, ylim[2] * 0.9, c("mixture", "fits", "appr"),
       col = c(1, 2, 4), lty = c(1, 2, 3), lwd = 1.5)
