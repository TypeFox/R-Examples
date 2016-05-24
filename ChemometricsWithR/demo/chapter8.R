## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("e1071")) {
  e1071.present <- FALSE
  cat("Package e1071 not available - some code may not run.\nInstall it by typing 'install.packages(\"e1071\")'")
} else {
  e1071.present <- TRUE
}
if (!require("nnet")) {
  nnet.present <- FALSE
  cat("Package nnet not available - some code may not run.\nInstall it by typing 'install.packages(\"nnet\")'")
} else {
  nnet.present <- TRUE
}

data(gasoline, package = "pls")
X <- gasoline$NIR[, 100*(1:4)]
Y <- gasoline$octane
odd <- seq(1, nrow(X), by = 2)
even <- seq(2, nrow(X), by = 2)
Xtr <- cbind(1, X[odd,])
Ytr <- Y[odd]
t(solve(crossprod(Xtr), t(Xtr)) %*% Ytr)

Xtr <- X[odd,]
Blm <- lm(Ytr ~ Xtr)
summary(Blm)

Xtr <- cbind(1, gasoline$NIR[odd,])
try(solve(crossprod(Xtr), t(Xtr)) %*% Ytr)

Blm <- ginv(Xtr) %*% Ytr

## PCR
X <- scale(gasoline$NIR, scale = FALSE,
           center = colMeans(gasoline$NIR[odd,]))
Xodd.svd <- svd(X[odd,])
Xodd.scores <- Xodd.svd$u %*% diag(Xodd.svd$d)
gasodd.pcr <- 
  lm(gasoline$octane[odd] ~ I(Xodd.scores[,1:5]) - 1)
gasodd.coefs <- coef(gasodd.pcr) %*% t(Xodd.svd$v[,1:5])
gasoline.pcr <- pcr(octane ~ ., data = gasoline, 
                    subset = odd, ncomp = 5)
all.equal(c(coef(gasoline.pcr)), c(gasodd.coefs))

wavelengths <- seq(900, 1700, by = 2)
plot(wavelengths, coef(gasoline.pcr), type = "l",
     xlab = "Wavelength (nm)", ylab = "Regression coefficient")

summary(gasoline.pcr)

RMSEP(gasoline.pcr, estimate = "train", intercept = FALSE)

RMSEP(gasoline.pcr, estimate = "train", comp = 4)

X <- gasoline$NIR[, 100*(1:4)]
Xtr <- X[odd,]
Blm <- lm(Ytr ~ Xtr)
rms(Ytr, fitted(Blm))

gasoline.pcr <- pcr(octane ~ ., data = gasoline, subset = odd,
                    validation = "LOO", ncomp = 10)
plot(gasoline.pcr, "validation", estimate = "CV")

RMSEP(gasoline.pcr, estimate = "CV")

sqrt(gasoline.pcr$validation$PRESS / nrow(Xtr))

par(pty = "s")
plot(gasoline.pcr, "prediction", ncomp = 4)
abline(0, 1, col = "gray")

gasoline.pcr.pred <- predict(gasoline.pcr, ncomp = 4,
                             newdata = gasoline[even,])
rms(gasoline$octane[even], gasoline.pcr.pred)

RMSEP(gasoline.pcr, ncomp = 4, newdata = gasoline[even,],
      intercept = FALSE)

## PLS
gasoline.pls <- plsr(octane ~ ., data = gasoline, 
                     subset = odd, ncomp = 5)
summary(gasoline.pls)

gasoline.pls <- plsr(octane ~ ., data = gasoline, subset = odd,
                     validation = "LOO", ncomp = 10)
plot(gasoline.pls, "validation", estimate = "CV")
opar <- par(pty = "s")
plot(gasoline.pls, "prediction", ncomp = 3)
abline(0, 1, col = "gray")
par(opar)

RMSEP(gasoline.pls, ncomp = 3, newdata = gasoline[even,],
      intercept = FALSE)

cor(gasoline.pls$loadings[,1:3])

cor(gasoline.pls$scores[,1:3])

plot(gasoline.pls, "loading", comps = 1:3, legendpos = "top", lwd = 1,
     lty = c(1, 2, 4), col = c(1, 2, 4), xlab = "variable nr")

plot(scores(gasoline.pls)[,1], Yscores(gasoline.pls)[,1],
     xlab = "X scores", ylab = "Y scores", main = "LV 1")
abline(h = 0, v = 0, col = "gray")
plot(scores(gasoline.pls)[,2], Yscores(gasoline.pls)[,2],
     xlab = "X scores", ylab = "Y scores", main = "LV 2")
abline(h = 0, v = 0, col = "gray")

## Ridge regression
gasoline.ridge <- 
  lm.ridge(octane ~ NIR, data = gasoline, subset = odd,
           lambda = seq(0.001, 0.1, by = 0.01))
select(gasoline.ridge)

if (e1071.present) {
  ## SVMs for regression
  gasoline.svm <- svm(octane ~ ., data = gasoline, 
                      subset = odd, cross = 10)
  
  plot(gasoline$octane[odd], predict(gasoline.svm),
       main = "Training set", xlab = "Octane number (true)", 
       ylab = "Octane number (predicted)")
  abline(0, 1)
  plot(gasoline$octane[even], 
       predict(gasoline.svm, new = gasoline[even,]),
       main = "Test set", xlab = "Octane number (true)", 
       ylab = "Octane number (predicted)")
  abline(0, 1)
  
  gasoline.svm <- svm(octane ~ ., data = gasoline, 
                      subset = odd, kernel = "linear")
  rms(gasoline$octane[even], 
      predict(gasoline.svm, new = gasoline[even,]))
}

if (nnet.present) {
  ## ANNs for regression
  X <- scale(gasoline$NIR, scale = FALSE,
             center = colMeans(gasoline$NIR[odd,]))
  Xodd.svd <- svd(X[odd,])
  Xodd.scores <- Xodd.svd$u %*% diag(Xodd.svd$d)
  Xeven.scores <- X[even,] %*% Xodd.svd$v
  
  set.seed(7)
  gas.nnet <- nnet(Xodd.scores[,1:5],
                   matrix(gasoline$octane[odd], ncol = 1),
                   size = 5, linout = TRUE)
  
  gas.nnet.pred <- predict(gas.nnet, Xeven.scores)
  rms(gas.nnet.pred, gasoline$octane[even])
}

data(wines, package = "ChemometricsWithRData")

## Classification by regression methods
wine.indices <- seq(1, 175, by = 25)
classvec2classmat(vintages[wine.indices])

## Error: these two lines appear twice in the book
C <- classvec2classmat(vintages[c(1:25, 61:85)])
X <- wines[c(1:25, 61:85), c(7, 13)]
wines.lm <- lm(C ~ X)
wines.lm.predict <- classmat2classvec(predict(wines.lm))
table(vintages[c(1:25, 61:85)], wines.lm.predict)

wines.lda <- lda(factor(vintages[c(1:25, 61:85)]) ~ X)
table(vintages[c(1:25, 61:85)], predict(wines.lda)$class)
