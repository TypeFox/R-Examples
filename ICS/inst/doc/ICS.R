### R code from vignette source 'ICS.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: ICS.Rnw:693-694
###################################################
options(digits = 4, prompt = "R> ", continue = "+  ")


###################################################
### code chunk number 2: ICS.Rnw:696-702
###################################################
library("ICS")
library("mvtnorm")
set.seed(2)
X <- rmvnorm(1000, c(0, 0, 1))
ics.X <- ics(X, stdKurt = FALSE)
Z <- ics.components(ics.X)


###################################################
### code chunk number 3: ICS.Rnw:705-706
###################################################
colMeans(X)


###################################################
### code chunk number 4: ICS.Rnw:710-711
###################################################
cov(X)


###################################################
### code chunk number 5: ICS.Rnw:715-716
###################################################
mean3(Z) - colMeans(Z)


###################################################
### code chunk number 6: ICS.Rnw:722-723
###################################################
(dim(X)[2] + 2) * (ics.X@gKurt - 1)


###################################################
### code chunk number 7: ICS.Rnw:745-754
###################################################
library("MASS")
library("ICS")
data("wood", package = "robustbase")
maha1.wood <- sqrt(mahalanobis(wood, colMeans(wood), cov(wood)))
set.seed(1)
covmve.wood <- cov.rob(wood)
maha2.wood <- sqrt(mahalanobis(wood, covmve.wood$center, covmve.wood$cov))
max.maha.wood <- max(c(maha1.wood, maha2.wood))
out.id <- ifelse(maha2.wood <= sqrt(qchisq(0.975, 6)), 0, 1)


###################################################
### code chunk number 8: mahalanobis (eval = FALSE)
###################################################
## par(mfrow = c(1, 2), las = 1)
## plot(maha1.wood, xlab = "index" ,ylab = "Mahalanobis distance", 
##   ylim = c(0, max.maha.wood), col = out.id + 1, pch = 15 * out.id + 1)
## abline(h = sqrt(qchisq(0.975, 6)))
## plot(maha2.wood, xlab = "index", ylab = "robust Mahalanobis distance", 
##   ylim = c(0, max.maha.wood), col = out.id + 1, pch = 15 * out.id + 1)
## abline(h = sqrt(qchisq(0.975, 6)))
## par(mfrow = c(1, 1))


###################################################
### code chunk number 9: ICS.Rnw:774-775
###################################################
par(mfrow = c(1, 2), las = 1)
plot(maha1.wood, xlab = "index" ,ylab = "Mahalanobis distance", 
  ylim = c(0, max.maha.wood), col = out.id + 1, pch = 15 * out.id + 1)
abline(h = sqrt(qchisq(0.975, 6)))
plot(maha2.wood, xlab = "index", ylab = "robust Mahalanobis distance", 
  ylim = c(0, max.maha.wood), col = out.id + 1, pch = 15 * out.id + 1)
abline(h = sqrt(qchisq(0.975, 6)))
par(mfrow = c(1, 1))


###################################################
### code chunk number 10: mahalanobis2 (eval = FALSE)
###################################################
## plot(maha1.wood, maha2.wood, xlab = "regular Mahalanobis distance",
##   ylab = "robust Mahalanobis distance", ylim = c(0, max.maha.wood),
##   xlim = c(0, max.maha.wood), col = out.id + 1, pch = 15 * out.id + 1,
##   las = 1)
## abline(0, 1)


###################################################
### code chunk number 11: ICS.Rnw:799-800
###################################################
plot(maha1.wood, maha2.wood, xlab = "regular Mahalanobis distance",
  ylab = "robust Mahalanobis distance", ylim = c(0, max.maha.wood),
  xlim = c(0, max.maha.wood), col = out.id + 1, pch = 15 * out.id + 1,
  las = 1)
abline(0, 1)


###################################################
### code chunk number 12: woodics (eval = FALSE)
###################################################
## library("ICSNP")
## my.HR.Mest <- function(X,...) HR.Mest(X,...)$scatter
## ics.default.wood <- ics(wood)
## ics.2.wood <- ics(wood, tM(wood)$V, tM(wood, 2)$V)
## ics.3.wood <- ics(wood, my.HR.Mest, HP1.shape)
## par(mfrow=c(1, 3), las = 1, mar = c(5, 4, 1, 1) + 0.1)
## plot(ics.components(ics.default.wood)[,6], xlab = "index", ylab = "IC 6",
##   sub = "ICS using cov and cov4", col = out.id + 1, pch = 15 * out.id + 1)
## plot(ics.components(ics.2.wood)[,6], xlab = "index", ylab = "IC 6",
##   sub = "ICS using tM(,1) and tM(,2)", col = out.id + 1,
##   pch = 15 * out.id + 1)
## plot(ics.components(ics.3.wood)[,6], xlab = "index", ylab = "IC 6",
##   sub = "ICS using HR.Mest and HP1.shape", col = out.id + 1,
##   pch = 15 * out.id + 1)
## par(mfrow = c(1, 1), las = 0)


###################################################
### code chunk number 13: ICS.Rnw:834-835
###################################################
library("ICSNP")
my.HR.Mest <- function(X,...) HR.Mest(X,...)$scatter
ics.default.wood <- ics(wood)
ics.2.wood <- ics(wood, tM(wood)$V, tM(wood, 2)$V)
ics.3.wood <- ics(wood, my.HR.Mest, HP1.shape)
par(mfrow=c(1, 3), las = 1, mar = c(5, 4, 1, 1) + 0.1)
plot(ics.components(ics.default.wood)[,6], xlab = "index", ylab = "IC 6",
  sub = "ICS using cov and cov4", col = out.id + 1, pch = 15 * out.id + 1)
plot(ics.components(ics.2.wood)[,6], xlab = "index", ylab = "IC 6",
  sub = "ICS using tM(,1) and tM(,2)", col = out.id + 1,
  pch = 15 * out.id + 1)
plot(ics.components(ics.3.wood)[,6], xlab = "index", ylab = "IC 6",
  sub = "ICS using HR.Mest and HP1.shape", col = out.id + 1,
  pch = 15 * out.id + 1)
par(mfrow = c(1, 1), las = 0)


###################################################
### code chunk number 14: irisplot (eval = FALSE)
###################################################
## library("ICS")
## library("MASS")
## data("iris")
## iris.ics <- ics(iris[,1:4])
## plot(iris.ics, col = as.numeric(iris[,5]))


###################################################
### code chunk number 15: ICS.Rnw:870-871
###################################################
library("ICS")
library("MASS")
data("iris")
iris.ics <- ics(iris[,1:4])
plot(iris.ics, col = as.numeric(iris[,5]))


###################################################
### code chunk number 16: irispairs (eval = FALSE)
###################################################
## pairs(princomp(iris[,1:4])$scores, col = as.numeric(iris[,5]))


###################################################
### code chunk number 17: ICS.Rnw:894-895
###################################################
pairs(princomp(iris[,1:4])$scores, col = as.numeric(iris[,5]))


###################################################
### code chunk number 18: irisics (eval = FALSE)
###################################################
## p <- dim(iris[, 1:4])[2]
## n <- dim(iris[, 1:4])[1]
## ngroup <- aggregate(iris$Species, list(iris$Species), length)$x
## colMeans.iris <- colMeans(iris[, 1:4])
## colMeans.iris.groups <- by(iris[, 1:4], iris$Species, colMeans)
## colMeans.iris.diffs <- sapply(colMeans.iris.groups,"-",
## colMeans.iris, simplify = FALSE)
## matrix.iris <- sapply(colMeans.iris.diffs, tcrossprod, simplify = FALSE)
## freq <- rep(ngroup, each = p^2)
## matrix.iris <- array(unlist(matrix.iris),
##   dim = c(p, p, nlevels(iris$Species)))
## cov.within <- rowSums(matrix.iris * freq, dims = 2)/n
## ics.iris.disc <- ics(iris[,1:4], cov(iris[,1:4]), cov.within)
## plot(ics.iris.disc, col = as.numeric(iris$Species))


###################################################
### code chunk number 19: ICS.Rnw:929-930
###################################################
p <- dim(iris[, 1:4])[2]
n <- dim(iris[, 1:4])[1]
ngroup <- aggregate(iris$Species, list(iris$Species), length)$x
colMeans.iris <- colMeans(iris[, 1:4])
colMeans.iris.groups <- by(iris[, 1:4], iris$Species, colMeans)
colMeans.iris.diffs <- sapply(colMeans.iris.groups,"-",
colMeans.iris, simplify = FALSE)
matrix.iris <- sapply(colMeans.iris.diffs, tcrossprod, simplify = FALSE)
freq <- rep(ngroup, each = p^2)
matrix.iris <- array(unlist(matrix.iris),
  dim = c(p, p, nlevels(iris$Species)))
cov.within <- rowSums(matrix.iris * freq, dims = 2)/n
ics.iris.disc <- ics(iris[,1:4], cov(iris[,1:4]), cov.within)
plot(ics.iris.disc, col = as.numeric(iris$Species))


###################################################
### code chunk number 20: iriskde (eval = FALSE)
###################################################
## iris.z <- ics.components(iris.ics)
## plot(density(iris.z[,4], bw = 0.15), las = 1,
## main = "Kernel Density of 4th component")
## rug(iris.z[1:50, 4], col = 1)
## rug(iris.z[51:100, 4], col = 2)
## rug(iris.z[101:150, 4], col = 3, ticksize = -0.03)


###################################################
### code chunk number 21: ICS.Rnw:955-956
###################################################
iris.z <- ics.components(iris.ics)
plot(density(iris.z[,4], bw = 0.15), las = 1,
main = "Kernel Density of 4th component")
rug(iris.z[1:50, 4], col = 1)
rug(iris.z[51:100, 4], col = 2)
rug(iris.z[101:150, 4], col = 3, ticksize = -0.03)


###################################################
### code chunk number 22: ICS.Rnw:971-982
###################################################
set.seed(4321)
train <- sample(1:150, 120)
lda.iris <- lda(Species ~ Sepal.Length + Sepal.Width + Petal.Length +
  Petal.Width, prior = c(1, 1, 1)/3, data = iris, subset = train)
table(iris[-train, 5], predict(lda.iris, iris[-train, ])$class)
ics.iris <- ics(as.matrix(iris[train, 1:4]))
iris.comp4 <- (ics.components(ics.iris))[,4]
lda.ics.iris <- lda(iris$Species[train] ~ iris.comp4, prior = c(1, 1, 1)/3)
iris.comp4.pred <- (as.matrix(iris[-train, 1:4]) %*% t(coef(ics.iris)))[,4]
table(iris[-train, 5], predict( lda.ics.iris, 
data.frame(iris.comp4 = iris.comp4.pred))$class)


###################################################
### code chunk number 23: ICS.Rnw:999-1001
###################################################
iris.centered <- sweep(iris[,1:4], 2, colMeans(iris[,1:4]), "-")
iris.paa <- ics(iris.centered, cov, covAxis, stdKurt = FALSE)


###################################################
### code chunk number 24: ICS.Rnw:1006-1009
###################################################
emp.align <- iris.paa@gKurt
mean(emp.align)
emp.align


###################################################
### code chunk number 25: screeplot (eval = FALSE)
###################################################
## screeplot(iris.paa, las = 1)
## abline(h = 1)


###################################################
### code chunk number 26: ICS.Rnw:1024-1025
###################################################
screeplot(iris.paa, las = 1)
abline(h = 1)


###################################################
### code chunk number 27: ICS.Rnw:1042-1047
###################################################
library("ICS")
library("pixmap")
fig1 <- read.pnm(system.file("pictures/cat.pgm", package = "ICS")[1])
fig2 <- read.pnm(system.file("pictures/road.pgm", package ="ICS")[1])
fig3 <- read.pnm(system.file("pictures/sheep.pgm", package = "ICS")[1])


###################################################
### code chunk number 28: ICS.Rnw:1054-1056
###################################################
p <- dim(fig1@grey)[2]
X <- cbind(as.vector(fig1@grey), as.vector(fig2@grey), as.vector(fig3@grey))


###################################################
### code chunk number 29: ICS.Rnw:1061-1065
###################################################
set.seed(4321)
A <- matrix(rnorm(9), ncol = 3)
X.mixed <- X %*% t(A)
ICA.fig <- ics(X.mixed, stdB="B")


###################################################
### code chunk number 30: pixmap (eval = FALSE)
###################################################
## par(mfrow = c(3, 3), omi = rep(0.1, 4), mai = rep(0.1, 4))
## plot(fig1)
## plot(fig2)
## plot(fig3)
## plot(pixmapGrey(X.mixed[,1], ncol = p))
## plot(pixmapGrey(X.mixed[,2], ncol = p))
## plot(pixmapGrey(X.mixed[,3], ncol = p))
## plot(pixmapGrey(ics.components(ICA.fig)[,1], ncol = p))
## plot(pixmapGrey(ics.components(ICA.fig)[,2], ncol = p))
## plot(pixmapGrey(ics.components(ICA.fig)[,3], ncol = p))


###################################################
### code chunk number 31: ICS.Rnw:1091-1094
###################################################
png(file = "ICS-pixmap.png", height=500, width=500)
par(mfrow = c(3, 3), omi = rep(0.1, 4), mai = rep(0.1, 4))
plot(fig1)
plot(fig2)
plot(fig3)
plot(pixmapGrey(X.mixed[,1], ncol = p))
plot(pixmapGrey(X.mixed[,2], ncol = p))
plot(pixmapGrey(X.mixed[,3], ncol = p))
plot(pixmapGrey(ics.components(ICA.fig)[,1], ncol = p))
plot(pixmapGrey(ics.components(ICA.fig)[,2], ncol = p))
plot(pixmapGrey(ics.components(ICA.fig)[,3], ncol = p))
dev.off()


###################################################
### code chunk number 32: ICS.Rnw:1115-1118
###################################################
library("ICS")
library("mvtnorm")
library("ICSNP")


###################################################
### code chunk number 33: ICS.Rnw:1125-1132
###################################################
set.seed(2000)
X <- rmvnorm(150, c(1, 2,-1))
A <- matrix(rnorm(9), ncol = 3)
b <- c(1, 1, 1)
X.trans <- sweep(X %*% t(A), 1, b, "+")
HL.estimator <- function(x){
wilcox.test(x, exact = TRUE, conf.int = TRUE)$estimate}


###################################################
### code chunk number 34: ICS.Rnw:1137-1140
###################################################
HLE.X <- apply(X, 2, HL.estimator)
as.vector(HLE.X %*% t(A) + b)
apply(X.trans, 2, HL.estimator)


###################################################
### code chunk number 35: ICS.Rnw:1151-1162
###################################################
ics.X <- ics(X, S1 = cov, S2 = tyler.shape)
HL.ics.Z1 <- apply(ics.components(ics.X), 2, HL.estimator)
HL.ics.X <- as.vector(HL.ics.Z1 %*% t(solve(coef(ics.X))))

ics.X.trans <- ics(X.trans, S1 = cov, S2 = tyler.shape)
HL.ics.Z2 <- apply(ics.components(ics.X.trans), 2, HL.estimator)
HL.ics.X.trans <- as.vector(HL.ics.Z2 %*% t(solve(coef(ics.X.trans))))

as.vector(HL.ics.X %*% t(A) +b)

HL.ics.X.trans


###################################################
### code chunk number 36: ICS.Rnw:1166-1170
###################################################
set.seed(1)
Y <- rmvt(60, diag(4), df = 6) + matrix(rep(c(0, 0.48), c(3*60, 60)),
  ncol = 4)
A2 <- matrix(rnorm(16), ncol = 4)


###################################################
### code chunk number 37: ICS.Rnw:1175-1177
###################################################
rank.ctest(Y, scores = "normal")
rank.ctest((Y %*% t(A2)), scores = "normal")


###################################################
### code chunk number 38: ICS.Rnw:1184-1190
###################################################
Z.Y <- as.matrix(ics.components(ics(Y,
  S1 = covOrigin, S2 = cov4, S2args = list(location = "Origin"))))
rank.ctest(Z.Y, scores = "normal")
Z.Y.trans <- as.matrix(ics.components(ics(Y %*% t(A2),
  S1 = covOrigin, S2 = cov4, S2args = list(location = "Origin"))))
rank.ctest(Z.Y.trans , scores = "normal")


