### R code from vignette source 'Multivariate_Extremes.Rnw'

###################################################
### code chunk number 1: rawdata
###################################################
options(show.signif.stars=FALSE)
library(evd); nn <- nrow(lossalae)
loss <- lossalae/1e+05; lts <- c(1e-04, 100)
plot(loss, log = "xy", xlim = lts, ylim = lts)


###################################################
### code chunk number 2: Multivariate_Extremes.Rnw:42-44
###################################################
ula <- apply(loss, 2, rank)/(nn + 1)
plot(ula)


###################################################
### code chunk number 3: Multivariate_Extremes.Rnw:49-50
###################################################
options(show.signif.stars=FALSE)
library(evd); nn <- nrow(lossalae)
loss <- lossalae/1e+05; lts <- c(1e-04, 100)
plot(loss, log = "xy", xlim = lts, ylim = lts)


###################################################
### code chunk number 4: asylogdfn
###################################################
abvevd(dep = 0.5, asy = c(1,1), model = "alog", plot = TRUE)
abvevd(dep = 0.5, asy = c(0.6,0.9), model = "alog", add = TRUE, lty = 2)
abvevd(dep = 0.5, asy = c(0.8,0.5), model = "alog", add = TRUE, lty = 3)


###################################################
### code chunk number 5: Multivariate_Extremes.Rnw:88-91
###################################################
abvevd(dep = -1/(-2), model = "neglog", plot = TRUE)
abvevd(dep = -1/(-1), model = "neglog", add = TRUE, lty = 2)
abvevd(dep = -1/(-0.5), model = "neglog", add = TRUE, lty = 3)


###################################################
### code chunk number 6: Multivariate_Extremes.Rnw:94-97
###################################################
abvevd(alpha = 1, beta = -0.2, model = "amix", plot = TRUE)
abvevd(alpha = 0.6, beta = 0.1, model = "amix", add = TRUE, lty = 2)
abvevd(alpha = 0.2, beta = 0.2, model = "amix", add = TRUE, lty = 3)


###################################################
### code chunk number 7: Multivariate_Extremes.Rnw:100-103
###################################################
abvevd(dep = 1/1.25, model = "hr", plot = TRUE)
abvevd(dep = 1/0.83, model = "hr", add = TRUE, lty = 2)
abvevd(dep = 1/0.5, model = "hr", add = TRUE, lty = 3)


###################################################
### code chunk number 8: Multivariate_Extremes.Rnw:108-109
###################################################
abvevd(dep = 0.5, asy = c(1,1), model = "alog", plot = TRUE)
abvevd(dep = 0.5, asy = c(0.6,0.9), model = "alog", add = TRUE, lty = 2)
abvevd(dep = 0.5, asy = c(0.8,0.5), model = "alog", add = TRUE, lty = 3)


###################################################
### code chunk number 9: Multivariate_Extremes.Rnw:121-129
###################################################
set.seed(131); cml <- loss[sample(nn),]
xx <- rep(1:50, each = 30); lts <- c(1e-04, 100)
cml <- cbind(tapply(cml[,1], xx, max), tapply(cml[,2], xx, max))
colnames(cml) <- colnames(loss)
plot(loss, log = "xy", xlim = lts, ylim = lts, col = "grey")
points(cml)
ecml <- -log(apply(cml,2,rank)/51)
plot(ecml)


###################################################
### code chunk number 10: nonpardfn
###################################################
pp <- "pickands"; cc <- "cfg"
abvnonpar(data = cml, epmar = TRUE, method = pp, plot = TRUE, lty = 3)
abvnonpar(data = cml, epmar = TRUE, method = pp, add = TRUE, madj = 1, lty = 2)
abvnonpar(data = cml, epmar = TRUE, method = pp, add = TRUE, madj = 2, lty = 4)
abvnonpar(data = cml, epmar = TRUE, method = cc, add = TRUE, lty = 1)


###################################################
### code chunk number 11: Multivariate_Extremes.Rnw:142-148
###################################################
m1 <- fbvevd(cml, asy1 = 1, model = "alog")
m2 <- fbvevd(cml, model = "log")
m3 <- fbvevd(cml, model = "bilog")
plot(m1, which = 4, nplty = 3)
plot(m2, which = 4, nplty = 3, lty = 2, add = TRUE)
plot(m3, which = 4, nplty = 3, lty = 4, add = TRUE)


###################################################
### code chunk number 12: Multivariate_Extremes.Rnw:153-154
###################################################
pp <- "pickands"; cc <- "cfg"
abvnonpar(data = cml, epmar = TRUE, method = pp, plot = TRUE, lty = 3)
abvnonpar(data = cml, epmar = TRUE, method = pp, add = TRUE, madj = 1, lty = 2)
abvnonpar(data = cml, epmar = TRUE, method = pp, add = TRUE, madj = 2, lty = 4)
abvnonpar(data = cml, epmar = TRUE, method = cc, add = TRUE, lty = 1)


###################################################
### code chunk number 13: Multivariate_Extremes.Rnw:164-167
###################################################
round(rbind(fitted(m2), std.errors(m2)), 3)
anova(m3, m2)
evind.test(cml, method = "score")


###################################################
### code chunk number 14: nonparqc
###################################################
lts <- c(0.01,100)
plot(loss, log = "xy", col = "grey", xlim = lts, ylim = lts)
points(cml); pp <- c(0.98,0.99,0.995)
qcbvnonpar(pp, data = cml, epmar = TRUE, mint = 30, add = TRUE)


###################################################
### code chunk number 15: Multivariate_Extremes.Rnw:185-186
###################################################
lts <- c(0.01,100)
plot(loss, log = "xy", col = "grey", xlim = lts, ylim = lts)
points(cml); pp <- c(0.98,0.99,0.995)
qcbvnonpar(pp, data = cml, epmar = TRUE, mint = 30, add = TRUE)


###################################################
### code chunk number 16: Multivariate_Extremes.Rnw:202-204 (eval = FALSE)
###################################################
## k0 <- bvtcplot(loss)$k0
## bvtcplot(loss, spectral = TRUE)


###################################################
### code chunk number 17: bvtc
###################################################
k0 <- bvtcplot(loss)$k0


###################################################
### code chunk number 18: Multivariate_Extremes.Rnw:213-214
###################################################
k0 <- bvtcplot(loss)$k0


###################################################
### code chunk number 19: Multivariate_Extremes.Rnw:224-228
###################################################
thresh <- apply(loss, 2, sort, decreasing = TRUE)[(k0+5)/2,]
mar1 <- fitted(fpot(loss[,1], thresh[1]))
mar2 <- fitted(fpot(loss[,2], thresh[2]))
rbind(mar1,mar2)


###################################################
### code chunk number 20: Multivariate_Extremes.Rnw:233-237
###################################################
m1 <- fbvpot(loss, thresh, model = "alog", asy1 = 1)
m2 <- fbvpot(loss, thresh, model = "bilog")
m3 <- fbvpot(loss, thresh, model = "bilog", likelihood = "poisson")
round(rbind(fitted(m2), std.errors(m2)), 3)


###################################################
### code chunk number 21: Multivariate_Extremes.Rnw:242-247
###################################################
abvnonpar(data = loss, method = "pot", k = k0, epmar = TRUE, 
  plot = TRUE, lty = 3)
plot(m1, which = 2, add = TRUE)
plot(m2, which = 2, add = TRUE, lty = 4)
plot(m3, which = 2, add = TRUE, lty = 2)


###################################################
### code chunk number 22: qcthresh
###################################################
lts <- c(1e-04, 100)
plot(loss, log = "xy", col = "grey", xlim = lts, ylim = lts)
plot(m1, which = 3, p = c(0.95,0.975,0.99), tlty = 0, add = TRUE)
abline(v=thresh[1], h=thresh[2])


###################################################
### code chunk number 23: Multivariate_Extremes.Rnw:261-262
###################################################
lts <- c(1e-04, 100)
plot(loss, log = "xy", col = "grey", xlim = lts, ylim = lts)
plot(m1, which = 3, p = c(0.95,0.975,0.99), tlty = 0, add = TRUE)
abline(v=thresh[1], h=thresh[2])


###################################################
### code chunk number 24: chiplot
###################################################
old <- par(mfrow = c(2,1))
chiplot(loss, ylim1 = c(-0.25,1), ylim2 = c(-0.25,1), nq = 200, 
  qlim = c(0.02,0.98), which = 1:2, spcases = TRUE)
par(old)


###################################################
### code chunk number 25: Multivariate_Extremes.Rnw:285-286
###################################################
old <- par(mfrow = c(2,1))
chiplot(loss, ylim1 = c(-0.25,1), ylim2 = c(-0.25,1), nq = 200, 
  qlim = c(0.02,0.98), which = 1:2, spcases = TRUE)
par(old)


###################################################
### code chunk number 26: etaplot
###################################################
fla <- apply(-1/log(ula), 1, min)
thresh <- quantile(fla, probs = c(0.025, 0.975))
tcplot(fla, thresh, nt = 100, pscale = TRUE, which = 2, vci = FALSE, 
  cilty = 2, type = "l", ylim = c(-0.2,1.2), ylab = "Tail Dependence")
abline(h = c(0,1))


###################################################
### code chunk number 27: Multivariate_Extremes.Rnw:304-307
###################################################
thresh <- quantile(fla, probs = 0.8)
m1 <- fpot(fla, thresh = thresh)
cat("Tail Dependence:", fitted(m1)["shape"], "\n")


###################################################
### code chunk number 28: Multivariate_Extremes.Rnw:310-312
###################################################
m2 <- fpot(fla, thresh = thresh, shape = 1)
anova(m1, m2, half = TRUE)


###################################################
### code chunk number 29: Multivariate_Extremes.Rnw:317-318
###################################################
fla <- apply(-1/log(ula), 1, min)
thresh <- quantile(fla, probs = c(0.025, 0.975))
tcplot(fla, thresh, nt = 100, pscale = TRUE, which = 2, vci = FALSE, 
  cilty = 2, type = "l", ylim = c(-0.2,1.2), ylab = "Tail Dependence")
abline(h = c(0,1))


