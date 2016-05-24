### R code from vignette source 'bdots.Rnw'

###################################################
### code chunk number 1: bdots.Rnw:34-36
###################################################
options(prompt = "R> ", continue = "+  ")
library("bdots")


###################################################
### code chunk number 2: bdots.Rnw:44-45 (eval = FALSE)
###################################################
## citation("bdots")


###################################################
### code chunk number 3: bdots.Rnw:58-73
###################################################
mini <- 0
peak <- 1
slope <- 0.002
cross <- 800

f <- function(time, mini, peak, slope, cross)
	mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini)))

plot(function(x) f(x, mini, peak, slope, cross), 0, 2000, xlab = "Time (ms)", ylab = "Proportion of Fixations")
lines(c(0, 500), c(0, 0), lty = 2)
lines(c(1000, 2000), c(1, 1), lty = 2)
points(800, 0.5, pch = 16)
lines(c(800, 800), c(0.2, 0.49), lty = 2)
lines(c(600, 1000), c(0.5 - 200 * slope, 0.5 + 200 * slope), lty = 2)
text(c(575, 925, 800, 1000), c(0, 1, 0.18, 0.92), c("B", "P", "C", "S"))


###################################################
### code chunk number 4: bdots.Rnw:86-118
###################################################
mu <- 800
ht <- 0.2
base1 <- 0
base2 <- 0.05
sig1 <- 175
sig2 <- 300

f <- function(time, mu, ht, base1, base2, sig1, sig2)
	(time < mu) * (exp(-1 * (time - mu) ^ 2 / (2 * sig1 ^ 2)) * (ht - base1) + base1) +
		(mu <= time) * (exp(-1 * (time - mu) ^ 2 / (2 * sig2 ^ 2)) * (ht - base2) + base2)

plot(function(x) f(x, mu, ht, base1, base2, sig1, sig2), 0, 2000, xlab = "Time (ms)", ylab = "Proportion of Fixations",
	xaxt = "n")
axis(1, at = c(0, 500, 800, 1000, 1500, 2000), labels = TRUE)
lines(c(0, 500), c(0, 0), lty = 2)
lines(c(1200, 2000), c(0.05, 0.05), lty = 2)
lines(c(500, 1100), c(0.20, 0.20), lty = 2)
lines(c(800, 800), c(0.02, 0.20), lty = 2)
lines(c(800, 800), c(-0.01, 0.01), lty = 2)
arrows(800, f(800 + 300 * 1.5, mu, ht, base1, base2, sig1, sig2),
	800 + 300 * 1.5, f(800 + 300 * 1.5, mu, ht, base1, base2, sig1, sig2),
	lty = 2, length = 0.12, code = 3)
arrows(800, f(800 - 175 * 1.5, mu, ht, base1, base2, sig1, sig2),
	800 - 175 * 1.5, f(800 - 175 * 1.5, mu, ht, base1, base2, sig1, sig2),
	lty = 2, length = 0.12, code = 3)
#text(c(625, 1075, 675, 1000, 800, 1200), c(0, 0.05, 0.07, 0.11, 0.015, 0.2),
#	c("Base 1", "Base 2", "Sig 1", "Sig 2", "Mu", "Height"))

text(c(560, 1150, 675, 1000, 800, 1150), c(0, 0.05, 0.07, 0.105, 0.014, 0.2),
	c(expression(B[1]), expression(B[2]),
		expression(sigma[1]), expression(sigma[2]),
		expression(mu), "H"))


###################################################
### code chunk number 5: bdots.Rnw:138-141
###################################################
data(ci)
names(ci)[1] <- "Group"
head(ci)


###################################################
### code chunk number 6: bdots.Rnw:154-155
###################################################
ci.1 <- subset(ci, ci$LookType == "Target")


###################################################
### code chunk number 7: bdots.Rnw:187-192
###################################################
ci.1.out.1 <- logistic.fit(ci.1, col = 4, rho.0 = 0.9, cor = TRUE,
	cores = 2)

ci.1.out.1$cor.1
ci.1.out.1$cor.2


###################################################
### code chunk number 8: bdots.Rnw:205-209
###################################################
ci.2 <- subset(ci, ci$LookType == "Cohort" | ci$LookType == "Unrelated")
ci.2$Curve <- ifelse(ci.2$LookType == "Cohort", 1, 2)
ci.2.out.1 <- doubleGauss.fit(ci.2, col = 4, rho.0 = 0.9, cor = TRUE, 
	cores = 1, diffs = TRUE)


###################################################
### code chunk number 9: bdots.Rnw:221-222 (eval = FALSE)
###################################################
## subs.plot(ci.1.out.1, "topleft")


###################################################
### code chunk number 10: bdots.Rnw:225-226
###################################################
ests.plot(ci.1.out.1)


###################################################
### code chunk number 11: bdots.Rnw:233-234 (eval = FALSE)
###################################################
## subs.plot(ci.2.out.1)


###################################################
### code chunk number 12: bdots.Rnw:237-238
###################################################
ests.plot(ci.2.out.1)


###################################################
### code chunk number 13: bdots.Rnw:259-265
###################################################
ci.2.out.1 <- doubleGauss.refit(ci.2.out.1 ,
	subj = c(13, 23),
	group = c(2, 2),
	curves = c(2, 2),
	params = list(c(650, 0.15, 150, 100, 0, 0.03),
		c(700, 0.10, 150, 100, 0, 0.01)))


###################################################
### code chunk number 14: bdots.Rnw:273-279
###################################################
refit.matrix <- matrix(NA, nrow = 2, ncol = 9)
refit.matrix[1,] <- c(13, 2, 2, 650, 0.15, 150, 100, 0, 0.03)
refit.matrix[2,] <- c(23, 2, 2, 700, 0.10, 150, 100, 0, 0.01)
refit.matrix
ci.2.out.1 <- doubleGauss.refit(ci.2.out.1, cor = FALSE,
	info.matrix = refit.matrix)


###################################################
### code chunk number 15: bdots.Rnw:282-287
###################################################
ci.2.out.1 <- doubleGauss.refit(ci.2.out.1 ,
	subj = c(13, 23),
	group = c(2, 2),
	curves = c(2, 2),
	cor = c(FALSE, FALSE))


###################################################
### code chunk number 16: bdots.Rnw:327-328
###################################################
ci.1.out.2 <- logistic.boot(ci.1.out.1, seed = 123, cores = 2)


###################################################
### code chunk number 17: bdots.Rnw:334-335
###################################################
logistic.boot(ci.1.out.1, seed = 123, cores = 2, test.params = TRUE)


###################################################
### code chunk number 18: bdots.Rnw:340-341
###################################################
ci.2.out.2 <- doubleGauss.boot(ci.2.out.1, seed = 123, cores = 2, time.test = 900)


###################################################
### code chunk number 19: bdots.Rnw:352-353
###################################################
replot(ci.1.out.2, bucket.lim = c(0, 1), main = "Example 1 Curve")


###################################################
### code chunk number 20: bdots.Rnw:357-359
###################################################
replot(ci.2.out.2, ylim = c(-0.01, 0.1), bucket.lim = c(0, 0.08),
	main = "Example 2 Curve")


