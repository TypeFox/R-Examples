## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("rrcov")) {
  rrcov.present <- FALSE
  cat("Package rrcov not available - some code may not run.\nInstall it by typing 'install.packages(\"rrcov\")'")
} else {
  rrcov.present <- TRUE
}
if (!require("ALS")) {
  ALS.present <- FALSE
  cat("Package ALS not available - some code may not run.\nInstall it by typing 'install.packages(\"ALS\")'")
} else {
  ALS.present <- TRUE
}
require(pls)

if (rrcov.present) {
  data(wines, package = "ChemometricsWithRData")
  ## Robust PCA
  X <- wines[vintages == "Grignolino",]
  X.sc <- scale(X)
  X.clPCA <- princomp(X.sc)
  X.robPCA <- princomp(X.sc, covmat = cov.mcd(X.sc))
  biplot(X.clPCA, main = "Classical PCA")
  biplot(X.robPCA, main = "MCD-based PCA")
  
  X.HubPCA5 <- PcaHubert(X.sc, k = 5)
  summary(X.HubPCA5)
  
  X.HubPCA <- PcaHubert(X.sc)
  summary(X.HubPCA)
  
  plot(X.HubPCA)
  
  ## OPLS
  data(gasoline, package = "pls")
  odd <- seq(1, length(gasoline$octane), by = 2)
  even <- seq(2, length(gasoline$octane), by = 2)
  gasoline$NIR <- scale(gasoline$NIR, scale = FALSE,
                        center = colMeans(gasoline$NIR[odd,]))
  Xtr <- gasoline$NIR[odd,]
  gasoline.pls <- plsr(octane ~ ., data = gasoline,
                       ncomp = 5, subset = odd,
                       validation = "LOO")
  ww <- gasoline.pls$loading.weights[,1]
  pp <- gasoline.pls$loadings[,1]
  w.ortho <- pp - crossprod(ww, pp)/crossprod(ww) * ww
  t.ortho <- Xtr %*% w.ortho
  p.ortho <- crossprod(Xtr, t.ortho) / c(crossprod(t.ortho))
  Xcorr <- Xtr - tcrossprod(t.ortho, p.ortho)
  gasoline.osc1 <- data.frame(octane = gasoline$octane[odd],
                              NIR = Xcorr)
  gasoline.opls1 <- plsr(octane ~ ., data = gasoline.osc1,
                         ncomp = 5, validation = "LOO")
  pp2 <- gasoline.opls1$loadings[,1]
  w.ortho2 <- pp2 - crossprod(ww, pp2)/crossprod(ww) * ww
  t.ortho2 <- Xcorr %*% w.ortho2
  p.ortho2 <- crossprod(Xcorr, t.ortho2) / c(crossprod(t.ortho2))
  Xcorr2 <- Xcorr - tcrossprod(t.ortho2, p.ortho2)
  gasoline.osc2 <- data.frame(octane = gasoline$octane[odd],
                              NIR = Xcorr2)
  gasoline.opls2 <- plsr(octane ~ ., data = gasoline.osc2,
                         ncomp = 5, validation = "LOO")
  plot(gasoline.pls, "validation", estimate = "CV",
       ylim = c(0.2, 1.5),
       main = "Gasoline training data (validation)")
  lines(0:5, c(RMSEP(gasoline.opls1, estimate = "CV"))$val,
        col = 2, lty = 2)
  lines(0:5, c(RMSEP(gasoline.opls2, estimate = "CV"))$val,
        col = 4, lty = 4)
  legend("topright", lty = c(1,2,4), col = c(1,2,4),
         legend = c("PLS", "OPLS: 1 OSC component", 
           "OPLS: 2 OSC components"))
  
  Xtst <- gasoline$NIR[even,]
  t.tst <- Xtst %*% w.ortho
  p.tst <- crossprod(Xtst, t.tst) / c(crossprod(t.tst))
  Xtst.osc1 <- Xtst - tcrossprod(t.tst, p.tst)
  gasoline.opls1.pred <- predict(gasoline.opls1,
                                 newdata = Xtst.osc1,
                                 ncomp = 2)
  t.tst2 <- Xtst.osc1 %*% w.ortho2
  p.tst2 <- crossprod(Xtst.osc1, t.tst2) / c(crossprod(t.tst2))
  Xtst.osc2 <- Xtst.osc1 - tcrossprod(t.tst2, p.tst2)
  gasoline.opls2.pred <- predict(gasoline.opls2,
                                 newdata = Xtst.osc2,
                                 ncomp = 1)
  RMSEP(gasoline.pls, newdata = gasoline[even,],
        ncomp = 3, intercept = FALSE)
  
  rms(gasoline$octane[even], gasoline.opls1.pred)
  
  rms(gasoline$octane[even], gasoline.opls2.pred)
  
  ## PCDA
  data(prostate, package = "ChemometricsWithRData")
  prost <- prostate[prostate.type != "bph",]
  prost.type <- factor(prostate.type[prostate.type != "bph"])
  
  prost.tcp <- tcrossprod(scale(prost))
  prost.svd <- svd(prost.tcp)
  prost.scores <- prost.svd$u %*% diag(sqrt(prost.svd$d))
  pairs(prost.scores[,1:5], pch = as.integer(prost.type),
        col = as.integer(prost.type),
        labels = paste("PC", 1:5))
  
  prost.pcda5 <- lda(prost.type ~ prost.scores[,1:5], CV = TRUE) ## INCORRECT
  table(prost.type, prost.pcda5$class)
  
  odd <- seq(1, length(prost.type), by = 2)
  even <- seq(2, length(prost.type), by = 2)
  prost.df <- data.frame(class = as.integer(prost.type),
                         msdata = I(prost))
  
  set.seed(7)
  prost.pcr <- pcr(class ~ msdata, ncomp = 16,
                   data = prost.df, subset = odd,
                   validation = "CV", scale = TRUE)
  validationplot(prost.pcr, estimate = "CV")
  
  prost.trn <- predict(prost.pcr)
  prost.trn.cl <- round(prost.trn[,1,])
  prost.trn.err <- apply(prost.trn.cl, 2, 
                         err.rate, prost.df$class[odd])
  plot(prost.trn.err, type = "l", col = 2,
       xlab = "# PCs", ylab = "Misclassif. rate")
  prost.tst <- predict(prost.pcr, newdata = prost.df[even,])
  prost.tst.cl <- round(prost.tst[,1,])
  prost.tst.err <- apply(prost.tst.cl, 2, 
                         err.rate, prost.df$class[even])
  lines(prost.tst.err, lty = 2)
  legend("topright", legend = c("LOO training set", "test set"),
         lty = c(1,2), col = 2:1)
  
  odd <- seq(1, length(prostate.type), by = 2)
  even <- seq(2, length(prostate.type), by = 2)
  prostate.clmat <- classvec2classmat(prostate.type)
  prostate.df <- data.frame(class = I(prostate.clmat),
                            msdata = I(prostate))
  
  set.seed(7)
  prostate.pcr <- pcr(class ~ msdata, ncomp = 16,
                      data = prostate.df, subset = odd,
                      validation = "CV", scale = TRUE)
  predictions.loo <-
    sapply(1:16, function(i, arr) classmat2classvec(arr[,,i]),
           prostate.pcr$validation$pred)
  loo.err <- apply(predictions.loo, 2, err.rate, 
                   prostate.type[odd])
  plot(loo.err, type = "l", main = "PCDA", col = 2, 
       ylim = c(.0, .47), xlab = "# PCs", ylab = "Misclassif. rate")
  
  prostate.pcrpred <- 
    predict(prostate.pcr, new = prostate.df[even,])
  predictions.pcrtest <-
    sapply(1:16, function(i, arr) classmat2classvec(arr[,,i]),
           prostate.pcrpred)
  lines(apply(predictions.pcrtest, 2, err.rate, 
              prostate.type[even]),
        type = "l", lty = 2)
  legend("bottomleft", legend = c("LOO training set", "test set"),
         lty = c(1,2), col = 2:1)
  
  table(prostate.type[even], predictions.pcrtest[,8])
  
  ## PLSDA
  set.seed(7)
  prostate.pls <- plsr(class ~ msdata, ncomp = 16,
                       data = prostate.df, subset = odd,
                       validation = "CV", scale = TRUE) 
  predictions.loo <-
    sapply(1:16, function(i, arr) classmat2classvec(arr[,,i]),
           prostate.pls$validation$pred)
  loo.err <- apply(predictions.loo, 2, err.rate, prostate.type[odd])
  prostate.plspred <- predict(prostate.pls, new = prostate.df[even,])
  predictions.plstest <-
    sapply(1:16, function(i, arr) classmat2classvec(arr[,,i]),
           prostate.plspred)
  table(prostate.type[even], predictions.plstest[,6])
  
  Xtst <- scale(prostate[even,],
                center = colMeans(prostate[odd,]),
                scale = apply(prostate[odd,], 2, sd))
  tst.scores <- Xtst %*% prostate.pls$projection
  prostate.ldapls <- lda(scores(prostate.pls)[,1:6], 
                         prostate.type[odd])
  table(prostate.type[even],
        predict(prostate.ldapls, new = tst.scores[,1:6])$class)
  
  plot(loo.err, type = "l", col = 2, ylim = c(.0, .47),
       xlab = "# LVs", ylab = "Misclassif. rate", main = "PLSDA")
  lines(apply(predictions.plstest, 2, err.rate, prostate.type[even]),
        type = "l", lty = 2)
  legend("bottomleft", legend = c("LOO training set", "test set"),
         lty = c(1,2), col = 2:1)
}
  
## Random data
nvar <- 2000
nobj <- 40
RandX <- matrix(rnorm(nobj*nvar), nrow = nobj)
RandY <- sample(c(0, 1), nobj, replace = TRUE)
Rand.pcr <- pcr(RandY ~ RandX, ncomp = 2)
Rand.ldapcr <- lda(RandY ~ scores(Rand.pcr), CV = TRUE)
table(RandY, Rand.ldapcr$class)

Rand.pls <- plsr(RandY ~ RandX, ncomp = 2)
Rand.ldapls <- lda(RandY ~ scores(Rand.pls), CV = TRUE)
table(RandY, Rand.ldapls$class)

## Calibration Transfer
data(shootout)
wl <- seq(600, 1898, by = 2)
indices <- which(wl >= 1100 & wl <= 1700)
nir.training1 <-
  data.frame(X = I(shootout$calibrate.1[,indices]), 
             y = shootout$calibrate.Y[,3])
mod1 <- plsr(y ~ X, data = nir.training1,
             ncomp = 5, validation = "LOO")
RMSEP(mod1, estimate = "CV")

RMSEP(mod1, estimate = "test", ncomp = 3, intercept = FALSE,
      newdata = data.frame(y = shootout$test.Y[,3], 
        X = I(shootout$test.1[,indices])))

RMSEP(mod1, estimate = "test", ncomp = 3, intercept = FALSE,
      newdata = data.frame(y = shootout$test.Y[,3], 
        X = I(shootout$test.2[,indices])))

nir.training2 <-
  data.frame(X = I(shootout$calibrate.2[,indices]),
             y = shootout$calibrate.Y[,3])
mod2 <- plsr(y ~ X, data = nir.training2,
             ncomp = 5, validation = "LOO")
plot(seq(1100, 1700, by = 2), coef(mod1, ncomp = 3), type = "l",
     xlab = "wavelength (nm)", ylab = "model coefficients",
     col = 1)
lines(seq(1100, 1700, by = 2), coef(mod2, ncomp = 3), col = 2)
legend("top", legend = c("set 1", "set 2"), bty = "n",
       col = 1:2, lty = 1)

recal.indices <- 1:5 * 10
F1 <- ginv(shootout$calibrate.2[recal.indices, indices]) %*%
  shootout$calibrate.1[recal.indices, indices]
RMSEP(mod1, estimate = "test", ncomp = 3, intercept = FALSE,
      newdata = data.frame(y = shootout$test.Y[,3], 
        X = I(shootout$test.2[,indices] %*% F1)))

contour(wl[indices], wl[indices], F1)

## MCR
data(bdata)
X <- bdata$d1
persp(X, phi = 20, theta = 34, expand = .5,
      xlab = "Time", ylab = "Wavelength")

X.efa <- efa(X, 3)
matplot(X.efa$forward, type = "l", col = 3:1,
        main = "Forward pass", ylab = "Singular values")
matplot(X.efa$backward, type = "l", col = 3:1,
        main = "Backward pass", ylab = "Singular values")

matplot(X.efa$pure.comp, type = "l", ylab = "", col = 3:1)

X.opa <- opa(X, 3)
matplot(X.opa$pure.comp, type = "l", col = 3:1,
        ylab = "response", xlab = "wavelength number")

X.mcr.efa <- mcr(X, X.efa$pure.comp, what = "col")
matplot(X.mcr.efa$C, type = "n", 
        main = "Concentration profiles",
        ylab = "Concentration")
matlines(X.efa$pure.comp, type = "l", lty = c(1,2,4), 
         col = "gray")
matlines(X.mcr.efa$C, type = "l", lty = c(1,2,4), col = 3:1)

matplot(t(X.mcr.efa$S), col = 3:1, type = "l", lty = c(1,2,4),
        main = "Pure spectra", ylab = "Intensity")

X.mcr.efa$rms

X.mcr.opa <- mcr(X, t(X.opa$pure.comp), what = "row")
X.mcr.opa$rms

if (ALS.present) {
  X.als.efa <- als(CList = list(X.efa$pure.comp),                          
                   PsiList = list(X), S = matrix(0, 73, 3),
                   nonnegS = TRUE, nonnegC = TRUE,
                   optS1st = TRUE, uniC = TRUE)
  X.als.opa <- als(CList = list(matrix(0, 40, 3)),
                   PsiList = list(X), S = X.opa$pure.comp,
                   nonnegS = TRUE, nonnegC = TRUE,
                   optS1st = FALSE, uniC = TRUE)
  Sefa <- sweep(X.als.efa$S, 2,
                apply(X.als.efa$S, 2, function(x) sqrt(sum(x^2))),
                FUN = "/")
  Sopa <- sweep(X.als.opa$S, 2,
                apply(X.als.opa$S, 2, function(x) sqrt(sum(x^2))),
                FUN = "/")
  Cefa <- sweep(X.als.efa$C[[1]], 2,
                apply(X.als.efa$C[[1]], 2, function(x) sqrt(sum(x^2))),
                FUN = "/")
  Copa <- sweep(X.als.opa$C[[1]], 2,
                apply(X.als.opa$C[[1]], 2, function(x) sqrt(sum(x^2))),
                FUN = "/")
  matplot(Sefa, type = "n", main = "Pure spectra (EFA)",
          ylab = "Intensity")
  abline(h = 0, col = "gray")
  matlines(Sefa, type = "l", lty = c(1,2,4), col = 3:1)
  
  matplot(Cefa, type = "n", 
          main = "Concentration profiles (EFA)",
          ylab = "Concentration")
  abline(h = 0, col = "gray")
  matlines(Cefa, lty = c(1,2,4), col = 3:1)
  
  matplot(Sopa, type = "n", main = "Pure spectra (OPA)",
          ylab = "Intensity")
  abline(h = 0, col = "gray")
  matlines(Sopa, type = "l", lty = c(1,2,4), col = 3:1)
  
  matplot(Copa, type = "n", 
          main = "Concentration profiles (OPA)",
          ylab = "Concentration")
  abline(h = 0, col = "gray")
  matlines(Copa, lty = c(1,2,4), col = 3:1)
  
  C0 <- matrix(0, 40, 3)
  X2.als.opa <- als(CList = list(C0, C0),
                    PsiList = list(bdata$d1, bdata$d2),
                    S = X.opa$pure.comp,               
                    nonnegS = TRUE, nonnegC = TRUE,
                    optS1st = FALSE, uniC = TRUE)
  
  resids2 <- X2.als.opa$S[,1:2] - cbind(c(bdata$sp1), c(bdata$sp2))
  apply(resids2, 2, function(x) sum(x^2))
  
  X.als.opa <- als(CList = list(matrix(0, 40, 3)),
                   PsiList = list(X), S = X.opa$pure.comp,
                   nonnegS = TRUE, nonnegC = TRUE,
                   optS1st = FALSE, uniC = TRUE)
  resids <- X.als.opa$S[,1:2] - cbind(c(bdata$sp1), c(bdata$sp2))
  apply(resids, 2, function(x) sum(x^2))
}
