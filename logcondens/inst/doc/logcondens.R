### R code from vignette source 'logcondens.Rnw'

###################################################
### code chunk number 1: init
###################################################
options(prompt = "R> ", continue = "+  ", width = 75, digits = 4)

if(require(logcondens) == FALSE){install.packages("logcondens"); library("logcondens")}

data("reliability")

set.seed(1)
n.sim <- 40
x.sim <- sort(rnorm(n.sim))


###################################################
### code chunk number 2: genSimul (eval = FALSE)
###################################################
## library("logcondens")
## set.seed(1)
## n.sim <- 40
## x.sim <-sort(rnorm(n.sim))


###################################################
### code chunk number 3: compSimul (eval = FALSE)
###################################################
## res <- logConDens(x.sim, smoothed = FALSE, print = FALSE)
## xs <- seq(-5, 5, by = 0.01)
## f.true <- dnorm(xs)
## par(las = 1, oma = c(0, 0, 3, 0), mar = c(3, 3.5, 0.5, 0.5), 
##     mfrow = c(1, 2))
## plot(res, which = "density", add.title = FALSE, legend.pos = "none")
## title(main = "Log-concave density estimation from i.i.d. data", 
##     outer = TRUE)
## mtext("Dashed vertical lines indicate knots of the log-density", 3, 
##     outer = TRUE)
## lines(xs, f.true, col = 4, lwd = 2)
## lines(density(x.sim), lwd = 2)
## legend("topleft", c("Kernel", "Estimated", "True"), lty = 1, lwd = 2,
##     col = c(1, 2, 4), bty = "n")
## plot(res, which = "log-density", add.title = FALSE, legend.pos = "none")
## lines(xs, log(f.true), col = 4, lwd = 2)


###################################################
### code chunk number 4: plotSimul
###################################################
res <- logConDens(x.sim, smoothed = FALSE, print = FALSE)
xs <- seq(-5, 5, by = 0.01)
f.true <- dnorm(xs)
par(las = 1, oma = c(0, 0, 3, 0), mar = c(3, 3.5, 0.5, 0.5), 
    mfrow = c(1, 2))
plot(res, which = "density", add.title = FALSE, legend.pos = "none")
title(main = "Log-concave density estimation from i.i.d. data", outer = TRUE)
mtext("Dashed vertical lines indicate knots of the log-density", 3, outer = TRUE)
lines(xs, f.true, col = 4, lwd = 2)
lines(density(x.sim), lwd = 2)
legend("topleft", c("Kernel", "Estimated", "True"), lty = 1, lwd = 2,
    col = c(1, 2, 4), bty = "n")
plot(res, which = "log-density", add.title = FALSE, legend.pos = "none")
lines(xs, log(f.true), col = 4, lwd = 2)


###################################################
### code chunk number 5: dispSummary1
###################################################
res <- logConDens(x.sim, smoothed = TRUE, print = FALSE)
summary(res)


###################################################
### code chunk number 6: dispSummary2
###################################################
evaluateLogConDens(xs = -1, res, which = 1:5)
quantilesLogConDens(ps = 0.5, res)


###################################################
### code chunk number 7: H1 (eval = FALSE)
###################################################
## n <- res$n
## xn <- res$xn
## ss <- sort(unique(c(x.sim, seq(min(x.sim), max(x.sim), length = 200))))
## H1 <- intF(ss, res)
## H2 <- intECDF(ss, xn)
## Ht <- H1 - H2
## ED <- 0 : (n + 1) / n
## 
## par(mfrow = c(2, 1), mar = c(0, 5, 2, 1), mex = 1, las = 1)
## plot(x.sim, res$Fhat, type = 'n', xaxt = 'n', ylab = 'Distribution 
##     functions')
## rug(x.sim); lines(x.sim, res$Fhat, col = 2, lwd = 1.5)
## lines(c(min(xn) - 10, xn, max(xn) + 10), ED, type = 's', lwd = 1.5)
## abline(v = res$knots, lty = 3); par(mar = c(4.5, 5, 0, 1))
## legend(-2.1, 1, c("empirical distribution function", expression("CDF based 
##     on "*hat(f)[m])), lty = 1, col = 1:2, lwd = 2, bty = "n")
## 
## plot(ss, Ht, type = 'n', xlab = "generated random sample x", 
##     ylab = "process H(t)", yaxt = "n")
## lines(ss, Ht, col = 2, lwd = 1.5)
## ax <- -c(0.02, 0.01, 0); axis(2, at = ax , labels = ax, cex.axis = 0.8)
## rug(x.sim); abline(v = res$knots, lty = 3); abline(h = 0, lty = 1)


###################################################
### code chunk number 8: plot2
###################################################
n <- res$n
xn <- res$xn
ss <- sort(unique(c(x.sim, seq(min(x.sim), max(x.sim), length = 200))))
H1 <- intF(ss, res)
H2 <- intECDF(ss, xn)
Ht <- H1 - H2
ED <- 0 : (n + 1) / n

par(mfrow = c(2, 1), mar = c(0, 5, 2, 1), mex = 1, las = 1)
plot(x.sim, res$Fhat, type = 'n', xaxt = 'n', ylab = 'Distribution functions')
rug(x.sim); lines(x.sim, res$Fhat, col = 2, lwd = 1.5)
lines(c(min(xn) - 10, xn, max(xn) + 10), ED, type = 's', lwd = 1.5)
abline(v = res$knots, lty = 3); par(mar = c(4.5, 5, 0, 1))
legend(-2.1, 1, c("empirical distribution function", expression("CDF based on "*hat(f)[m])), lty = 1, col = 1:2, lwd = 2, bty = "n")

plot(ss, Ht, type = 'n', xlab = "t", ylab = "process H(t)", yaxt = "n")
lines(ss, Ht, col = 2, lwd = 1.5)
ax <- -c(0.02, 0.01, 0); axis(2, at = ax , labels = ax, cex.axis = 0.8)
rug(x.sim); abline(v = res$knots, lty = 3); abline(h = 0, lty = 1)


###################################################
### code chunk number 9: disp_rel1 (eval = FALSE)
###################################################
## x.rel <- sort(reliability)
## n <- length(x.rel)
## mu <- mean(x.rel); sig <- sd(x.rel)
## xs <- seq(1350, 1950, length.out = 500)
## res <- logConDens(x.rel, smoothed = TRUE, print = FALSE, xs = xs)
## f.smoothed <- res$f.smoothed
## xs2 <- xs[(xs >= min(x.rel)) & (xs <= max(x.rel))]
## f <- rep(NA, length(xs2))
## for (i in 1:length(xs2)){f[i] <- evaluateLogConDens(xs2[i], 
##     res)[, "density"]}
## h <- sig / sqrt(n)
## f.kernel <- rep(NA, length(xs))
## for (i in 1:length(xs)){f.kernel[i] <- mean(dnorm(xs[i], mean = 
##     x.rel, sd = h))}
## f.normal <- dnorm(xs, mean = mu, sd = sig)
## par(las = 1, mar = c(3, 3.5, 0.5, 0.5))
## plot(0, 0, type = 'n', xlim = c(1390, 1900), ylim = 
##     c(0, 6.5 * 10^-3), ylab = "")
## rug(x.rel)
## lines(xs, f.normal, col = 3)
## lines(xs, f.kernel, col = 4)
## lines(xs, f.smoothed, lwd = 4, col = 5)
## lines(xs2, f, col = 2)
## segments(c(-1300, max(x.rel)), c(0, 0), c(min(x.rel), 2000), 
##     c(0, 0), col = 2)
## legend("topleft", c(expression("log-concave "*hat(f)[n]),
##     expression("normal "*hat(f)[nor]), expression("kernel "*hat(f)[ker]),
##     expression("log-concave smoothed "*hat(f)[n]*"*")),
##     lty = 1, lwd = 3, col = 2:5, bty = "n")
## segments(res$knots, 0, res$knots, 0.002, lty = 2)


###################################################
### code chunk number 10: plot1
###################################################
x.rel <- sort(reliability)
n <- length(x.rel)
mu <- mean(x.rel); sig <- sd(x.rel)
xs <- seq(1350, 1950, length.out = 500)
res <- logConDens(x.rel, smoothed = TRUE, print = FALSE, xs = xs)
f.smoothed <- res$f.smoothed
xs2 <- xs[(xs >= min(x.rel)) & (xs <= max(x.rel))]
f <- rep(NA, length(xs2))
for (i in 1:length(xs2)){f[i] <- evaluateLogConDens(xs2[i], res)[, "density"]}
h <- sig / sqrt(n)
f.kernel <- rep(NA, length(xs))
for (i in 1:length(xs)){f.kernel[i] <- mean(dnorm(xs[i], mean = x.rel, sd = h))}
f.normal <- dnorm(xs, mean = mu, sd = sig)
par(las = 1, mar = c(3, 3.5, 0.5, 0.5))
plot(0, 0, type = 'n', xlim = c(1390, 1900), ylim = c(0, 6.5 * 10^-3), ylab = "")
rug(x.rel)
lines(xs, f.normal, col = 3)
lines(xs, f.kernel, col = 4)
lines(xs, f.smoothed, lwd = 4, col = 5)
lines(xs2, f, col = 2)
segments(c(-1300, max(x.rel)), c(0, 0), c(min(x.rel), 2000), c(0, 0), col = 2)
legend("topleft", c(expression("log-concave "*hat(f)[n]),
    expression("normal "*hat(f)[nor]), expression("kernel "*hat(f)[ker]),
    expression("log-concave smoothed "*hat(f)[n]*"*")),
    lty = 1, lwd = 3, col = 2:5, bty = "n")
segments(res$knots, 0, res$knots, 0.002, lty = 2)


###################################################
### code chunk number 11: setdig
###################################################
options(digits = 7)
Msimul <- 1000


###################################################
### code chunk number 12: simulrel1
###################################################
set.seed(1977)
rel_samples <- rlogcon(n = 20, x0 = x.rel)


###################################################
### code chunk number 13: simulrel2
###################################################
sort(rel_samples$X)


###################################################
### code chunk number 14: simulrel3
###################################################
sort(rel_samples$X_star)


###################################################
### code chunk number 15: setdig
###################################################
options(digits = 4)


###################################################
### code chunk number 16: disp2sample
###################################################
set.seed(1)
n1 <- 20
n2 <- 25
x <- sort(rgamma(n1, 2, 1))
y <- sort(rgamma(n2, 2, 1) + 0.5)
twosample <- logconTwoSample(x, y, M = 5, display = FALSE)
twosample$p.value
twosample$test.stat.orig
twosample$test.stats[1:5, ]


###################################################
### code chunk number 17: sim1
###################################################

#path <- "C:/rufibach/research/p01_math/p17_logcondens_package/paper/"
path <- ""
path.res  <- paste(path, "results/", sep = "")
par(mfrow = c(3, 1), oma = rep(0, 4), mar = c(4.5, 4.5, 2, 1))

## ------------------------------
## setting 1
## ------------------------------
alpha <- 0.05
mus <- c(0, 0.5, 1, 1.5, 2)
n1 <- 20
n2 <- 25
M0 <- 999
setting <- 1
mu1 <- 0

# number of decisicions on H_1 (which is true)
abs <- matrix(NA, nrow = length(mus), ncol = 10)
rel <- matrix(NA, nrow = length(mus), ncol = 7)
for (i in 1:length(mus)){
    
    mu2 <- mus[i]

    name <- paste(path.res, "pvals_setting=", setting, "_n1=", n1, "_n2=", n2, "_mu1=", mu1, "_mu2=", mu2, ".txt", sep='')
    p.vals <- read.table(file = name, sep = ",", header = TRUE)
    abs[i, 1:7] <- c(n1, n2, mu1, mu2, apply(p.vals <= alpha, 2, sum))
    abs[i, 8:10] <- apply(is.na(p.vals) == FALSE, 2, sum, na.rm = TRUE)

    rel[i, 1:7] <- c(n1, n2, mu1, mu2, apply(p.vals <= alpha, 2, sum) / Msimul)
}
 
abs <- as.data.frame(abs)
dimnames(abs)[[2]] <- c("n1", "n2", "mu1", "mu2", "H1 logcon", "H1 logcon smooth", "H1 K-S", "#tests Logcon", "#tests logcon smooth", "#tests K-S")
rel <- as.data.frame(rel)
dimnames(rel)[[2]] <- dimnames(abs)[[2]][1:7]
absS1 <- abs
relS1 <- rel
 

plot(0, 0, type = 'n', xlim = c(-0.1, 2.1), ylim = c(0, 1), xlab = expression("location difference "*mu), ylab = "proportion of rejected null hypothesis")
title(expression("Results for Setting 1:  N(0, 1) vs. N("*mu*", 1)"))
move <- 0.015
points(rel[, 4] - move, rel[, 5], col = 1, pch = 3, lwd = 2)      ## logcon test
points(rel[, 4], rel[, 6], col = 3, pch = 5, lwd = 2)             ## logcon smooth test
points(rel[, 4], rel[, 7], col = 2, pch = 4, lwd = 2)             ## K-S test
abline(h = c(1, alpha, 0), lty = 2)
legend(-0.1, 0.9, c("based on log-concave CDF", "Kolmogorov-Smirnov"), pch = c(3, 4), col = 1:2, pt.lwd = 2, title = "Type of test:", bty = "n")


## ------------------------------
## setting 2
## ------------------------------
alpha <- 0.05
n1 <- 20
n2 <- 25
mus <- c(0, 0.5, 1, 1.5, 2)
setting <- 2
mu1 <- 0

# number of decisicions on H_1 (which is true)
abs <- matrix(NA, nrow = length(mus), ncol = 10)
rel <- matrix(NA, nrow = length(mus), ncol = 7)
for (i in 1:length(mus)){
    
    mu2 <- mus[i]

    name <- paste(path.res, "pvals_setting=", setting, "_n1=", n1, "_n2=", n2, "_mu1=", mu1, "_mu2=", mu2, ".txt", sep='')
    p.vals <- read.table(file = name, sep = ",", header = TRUE)
    abs[i, 1:7] <- c(n1, n2, mu1, mu2, apply(p.vals <= alpha, 2, sum))
    abs[i, 8:10] <- apply(is.na(p.vals) == FALSE, 2, sum, na.rm = TRUE)

    rel[i, 1:7] <- c(n1, n2, mu1, mu2, apply(p.vals <= alpha, 2, sum) / Msimul)
}
 
dimnames(abs)[[2]] <- dimnames(absS1)[[2]]
dimnames(rel)[[2]] <- dimnames(relS1)[[2]]
absS2 <- abs
relS2 <- rel
 

plot(0, 0, type = 'n', xlim = c(-0.1, 2.1), ylim = c(0, 1), xlab = expression("location difference "*mu), ylab = "proportion of rejected null hypothesis")
title(expression("Results for Setting 2:  Gam(2, 1) vs. Gam(2, 1) + "*mu))
move <- 0.015
points(rel[, 4] - move, rel[, 5], col = 1, pch = 3, lwd = 2)      ## logcon test
points(rel[, 4], rel[, 6], col = 3, pch = 5, lwd = 2)             ## logcon smooth test
points(rel[, 4], rel[, 7], col = 2, pch = 4, lwd = 2)             ## K-S test
abline(h = c(1, alpha, 0), lty = 2)


## ------------------------------
## setting 3
## ------------------------------
alpha <- 0.05
n1 <- 20
n2 <- 25
taus <- c(1, 1.5, 2, 2.5, 3)
setting <- 3

# number of decisicions on H_1 (which is true)
abs <- matrix(NA, nrow = length(mus), ncol = 9)
rel <- matrix(NA, nrow = length(mus), ncol = 6)
for (i in 1:length(taus)){
    
    tau <- taus[i]
    
    name <- paste(path.res, "pvals_setting=", setting, "_n1=", n1, "_n2=", n2, "_tau=", tau, ".txt", sep='')   
    p.vals <- read.table(file = name, sep = ",", header = TRUE)
    abs[i, 1:6] <- c(n1, n2, tau, apply(p.vals <= alpha, 2, sum))
    abs[i, 7:9] <- apply(is.na(p.vals) == FALSE, 2, sum, na.rm = TRUE)

    rel[i, 1:6] <- c(n1, n2, tau, apply(p.vals <= alpha, 2, sum) / Msimul)
}
 
abs <- as.data.frame(abs)
dimnames(abs)[[2]] <- c("n1", "n2", "tau", "H1 logcon", "H1 logcon smooth", "H1 K-S", "#tests Logcon", "#tests logcon smooth", "#tests K-S")
rel <- as.data.frame(rel)
dimnames(rel)[[2]] <- dimnames(abs)[[2]][1:6]
absS3 <- abs
relS3 <- rel


plot(0, 0, type = 'n', xlim = c(0.9, 3.1), ylim = c(0, 1), xlab = expression("alternative shape parameter "*tau), ylab = "proportion of rejected null hypothesis")
title(expression("Results for Setting 3:  Gam(2, 1) vs. Gam("*tau*", 1)"))
move <- 0.015
points(rel[, 3] - move, rel[, 4], col = 1, pch = 3, lwd = 2)      ## logcon test
points(rel[, 3], rel[, 5], col = 3, pch = 5, lwd = 2)             ## logcon smooth test
points(rel[, 3], rel[, 6], col = 2, pch = 4, lwd = 2)             ## K-S test
abline(h = c(1, alpha, 0), lty = 2)


###################################################
### code chunk number 18: plot2sample
###################################################
par(mar = c(4.5, 4, 1, 1), las = 1)
set.seed(111)
n1 <- 20
n2 <- 25
x <- sort(rgamma(n1, 2, 1))
y <- sort(rgamma(n2, 2, 1) + 0.5)
res1 <- activeSetLogCon(x)
res2 <- activeSetLogCon(y)

grid <- unique(sort(c(seq(0, max(x, y), length.out = 500), x, y)))
Fhat <- matrix(NA, ncol = 2, nrow = length(grid))
for (i in 1:nrow(Fhat)){
    Fhat[i, 1] <- evaluateLogConDens(grid[i], res1, which = 3)[4]
    Fhat[i, 2] <- evaluateLogConDens(grid[i], res2, which = 3)[4]
}

plot(c(0, x), 0:n1 / n1, type = 's', xlim = c(range(c(x, y)) + c(0, 1)), ylim = c(0, 1), main = "", xlab = "data", ylab = "distribution functions", lwd = 2)
lines(c(0, y), 0:n2 / n2, type = 's', col = 2, lwd = 2)
lines(grid, Fhat[, 1], lwd = 2); rug(x, col = 1, lwd = 2)
lines(grid, Fhat[, 2], col = 2, lwd = 2); rug(y, col = 2, lwd = 2)
segments(c(0, max(x)), c(0, 1), c(min(x), 20), c(0, 1), col = 1, lwd = 2)
segments(c(0, max(y)), c(0, 1), c(min(y), 20), c(0, 1), col = 2, lwd = 2)
md <- maxDiffCDF(res1, res2, which = c("MLE", "smooth")[1])
stat <- md$test.stat[1]
loc <- md$location[1]
abline(v = loc, lty = 4)
abline(h = c(evaluateLogConDens(loc, res1, which = 3)[, "CDF"], evaluateLogConDens(loc, res2, which = 3)[, "CDF"]), lty = 4)

legend(4, 0.6, c(expression(Gamma*"(2, 1)"), expression(Gamma*"(2, 1) + 0.5")), lty = 1, col = 1:2, bty = "n", lwd = 3)


