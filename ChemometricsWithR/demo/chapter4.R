## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("fastICA")) {
  fastICA.present <- FALSE
  cat("Package fastICA not available - some code may not run.\nInstall it by typing 'install.packages(\"fastICA\")'")
} else {
  fastICA.present <- TRUE
}


data(wines, package = "ChemometricsWithRData")
wines.sc <- scale(wines)
## PCA machinery
wines.svd <- svd(wines.sc)
wines.scores <- wines.svd$u %*% diag(wines.svd$d)
wines.loadings <- wines.svd$v
wines.vars <- wines.svd$d^2 / (nrow(wines) - 1)
wines.totalvar <- sum(wines.vars)
wines.relvars <- wines.vars / wines.totalvar
variances <- 100 * round(wines.relvars, digits = 3)
variances[1:5]

plot(wines.scores[,1:2], type = "n",
     xlab = paste("PC 1 (", variances[1], "%)", sep = ""),
     ylab = paste("PC 2 (", variances[2], "%)", sep = ""))
abline(h = 0, v = 0, col = "gray")
points(wines.scores[,1:2], pch = wine.classes, col = wine.classes)
plot(wines.loadings[,1] * 1.2, wines.loadings[,2], type = "n",
     xlab = paste("PC 1 (", variances[1], "%)", sep = ""),
     ylab = paste("PC 2 (", variances[2], "%)", sep = ""))
arrows(0, 0, wines.loadings[,1], wines.loadings[,2],
       length = .15, angle = 20, col = "blue")
text(wines.loadings[,1:2], labels = colnames(wines))
 
data(prostate, package = "ChemometricsWithRData")
cat("Speed test for PCA: direct SVD first.\n")
system.time({
  prost.svd <- svd(prostate)
  prost.scores <- prost.svd$u %*% diag(prost.svd$d)
  prost.variances <- prost.svd$d^2 / (nrow(prostate) - 1)
  prost.loadings <- prost.svd$v
})

cat("PCA speed test part 2: SVD on a crossproduct matrix.\n")
system.time({
  prost.tcp <- tcrossprod(prostate)
  prost.svd <- svd(prost.tcp)
  prost.scores <- prost.svd$u %*% diag(sqrt(prost.svd$d))
  prost.variances <- prost.svd$d / (nrow(prostate) - 1)
  prost.loadings <- solve(prost.scores, prostate)
})

barplot(wines.vars[1:10], main = "Variances",
        names.arg = paste("PC", 1:10))
barplot(log(wines.vars[1:10]), main = "log(Variances)",
        names.arg = paste("PC", 1:10))
barplot(wines.relvars[1:10], main = "Relative variances",
        names.arg = paste("PC", 1:10))
barplot(cumsum(100 * wines.relvars[1:10]), 
        main = "Cumulative variances (%)",
        names.arg = paste("PC", 1:10), ylim = c(0, 100))

llambdas <- log(wines.vars)
CIwidth <- qnorm(.975) * sqrt(2 / (nrow(wines) - 1))
CIs <- cbind(exp(llambdas - CIwidth),
             wines.vars,
             exp(llambdas + CIwidth))
colnames(CIs) <- c("CI 0.025", "  Estimate", "  CI 0.975")
CIs[1:5,]

small.ones <- wines.vars[11:13]
n <- nrow(wines)
nsmall <- length(small.ones)
geo.mean <- prod(small.ones)^{1/nsmall}
mychisq <- (n - 1) * nsmall * log(mean(small.ones) / geo.mean)
ndf <- (nsmall + 2) * (nsmall - 1) / 2
1 - pchisq(mychisq, ndf)

X1 <- scale(wines[1:88,])
X1.svd <- svd(X1)
X1.pca <- list(scores = X1.svd$u %*% diag(X1.svd$d),
               loadings = X1.svd$v)
X2 <- scale(wines[89:177,],
            center = attr(X1, "scaled:center"),
            scale = attr(X1, "scaled:scale"))
X2.scores <- X2 %*% X1.pca$loadings
plot(rbind(X1.pca$scores, X2.scores),
     pch = rep(c(1,2), c(88, 89)), 
     col = rep(c(1,2), c(88, 89)),
     xlab = "PC 1", ylab = "PC 2")

odd <- seq(1, nrow(wines), by = 2)
even <- seq(2, nrow(wines), by = 2)
X1 <- scale(wines[odd,])
X2 <- scale(wines[even,],
            center = attr(X1, "scaled:center"),
            scale = attr(X1, "scaled:scale"))
## The next commands are not shown again in the book
X1.svd <- svd(X1)
X1.pca <- list(scores = X1.svd$u %*% diag(X1.svd$d),
               loadings = X1.svd$v)
X2.scores <- X2 %*% X1.pca$loadings
plot(rbind(X1.pca$scores, X2.scores),
     pch = rep(c(1,2), c(88, 89)), 
     col = rep(c(1,2), c(88, 89)),
     xlab = "PC 1", ylab = "PC 2")

data(gasoline, package = "pls")
  
nir.prcomp <- prcomp(gasoline$NIR)
summary(nir.prcomp, digits = 2)

plot(nir.prcomp)

nir.loadings <- nir.prcomp$rotation[,1:4]
offset <- c(0, 0.09) # to create space for labels
plot(nir.loadings[,1:2], type = "l",
     xlim = range(nir.loadings[,1]) + offset,
     xlab = "PC 1 (72.6%)", ylab = "PC 2 (11.3%)")
points(nir.loadings[c(386, 396), 1:2])
text(nir.loadings[c(386, 396), 1:2], pos = 4,
     labels = paste(c(1670, 1690), "nm"))
offset <- c(-0.12, 0.12) # to create space for labels
plot(nir.loadings[,3:4], type = "l",
     xlim = range(nir.loadings[,3]) + offset,
     xlab = "PC 3 (6.95%)", ylab = "PC 4 (4.60%)")
points(nir.loadings[c(154, 370, 398), 3:4])
text(nir.loadings[c(154, 370), 3:4], pos = 4,
     labels = paste(c(1206, 1638), "nm"))
text(nir.loadings[398, 3:4, drop = FALSE],
     labels = "1694 nm", pos = 2)

biplot(nir.prcomp)

extremes <- c(15,41,45,57)
## Error in the book: one closing bracket too many
Xextr <- scale(gasoline$NIR, scale = FALSE)[extremes,]
wavelengths <- seq(900, 1700, by = 2)
matplot(wavelengths, t(Xextr),
        type = "l", xlab = "Wavelength (nm)",
        ylab = "Intensity (mean-scaled)", lty = c(1,1,2,2),
        col = c(1, 2, 1, 2))
legend("bottomleft", legend = paste("sample", extremes),
       lty = c(1,1,2,2), col = c(1,2,1,2), bty = "n")

wines.PCA <- PCA(scale(wines))
scoreplot(wines.PCA, pch = wine.classes, col = wine.classes)
loadingplot(wines.PCA, show.names = TRUE)

## Multidimensional scaling
wines.dist <- dist(scale(wines))
wines.cmdscale <- cmdscale(wines.dist)
plot(wines.cmdscale, 
     pch = wine.classes, col = wine.classes, 
     main = "Principal Coordinate Analysis",
     xlab = "Coord 1", ylab = "Coord 2")

wines.sammon <- sammon(wines.dist)
wines.sammon <- sammon(wines.dist, magic = .00003)
plot(wines.sammon$points, main = "Sammon mapping",
     col = wine.classes, pch = wine.classes, 
     xlab = "Coord 1", ylab = "Coord 2")

wines.isoMDS <- isoMDS(wines.dist)
plot(wines.isoMDS$points, main = "Non-metric MDS",
     col = wine.classes, pch = wine.classes, 
     xlab = "Coord 1", ylab = "Coord 2")

if (fastICA.present) {
  ## ICA and Projection Pursuit
  ## Components depend on the setting of the random seed
  set.seed(7)
  wines.ica <- fastICA(wines.sc, 3)
  pairs(wines.ica$S, main = "ICA components",
        col = wine.classes, pch = wine.classes)
  
  ## the order of the components is arbitrary! Compared to the book, ICs
  ## 1 and 3 are swapped here. Just showing three out of five components
  ## can be very dangerous...
  set.seed(13)
  wines.ica5 <- fastICA(wines.sc, 5)
  pairs(wines.ica5$S[,1:3], 
        main = "ICA components (3 out of 5)",
        col = wine.classes, pch = wine.classes)
}

## Factor analysis
wines.fa <- factanal(wines.sc, 3, scores = "regression")
wines.fa
pairs(wines.fa$scores, main = "FA components",
      col = wine.classes, pch = wine.classes)
