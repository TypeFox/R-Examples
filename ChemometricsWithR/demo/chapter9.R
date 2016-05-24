## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("boot")) {
  boot.present <- FALSE
  cat("Package boot not available - some code may not run.\nInstall it by typing 'install.packages(\"boot\")'")
} else {
  boot.present <- TRUE
}
if (!require("ipred")) {
  ipred.present <- FALSE
  cat("Package ipred not available - some code may not run.\nInstall it by typing 'install.packages(\"ipred\")'")
} else {
  ipred.present <- TRUE
}
if (!require("randomForest")) {
  randomForest.present <- FALSE
  cat("Package randomForest not available - some code may not run.\nInstall it by typing 'install.packages(\"randomForest\")'")
} else {
  randomForest.present <- TRUE
}
if (!require("ada")) {
  ada.present <- FALSE
  cat("Package ada not available - some code may not run.\nInstall it by typing 'install.packages(\"ada\")'")
} else {
  ada.present <- TRUE
}

data(gasoline, package = "pls")
gasoline.mscpcr <- pcr(octane ~ msc(NIR), data = gasoline, 
                       ncomp = 4)
gasoline.mscpcr.cv <- crossval(gasoline.mscpcr, loo = TRUE)
RMSEP(gasoline.mscpcr.cv, estimate = "CV")

data(wines, package = "ChemometricsWithRData")
X <- wines[vintages != "Barolo", ]
vint <- factor(vintages[vintages != "Barolo"])
kvalues <- 1:12
ktabs <- lapply(kvalues,
                function(i) {
                  kpred <- knn.cv(X, vint, k = i)
                  table(vint, kpred)
                })
TPrates <- sapply(ktabs, function(x) x[1,1]/sum(x[,1]))
FPrates <- sapply(ktabs, function(x) 1 - x[2,2]/sum(x[2,]))
plot(FPrates, TPrates, type = "b",
     xlim = c(.15, .45), ylim = c(.5, .75),
     xlab = "FP rate", ylab = "TP rate")
text(FPrates, TPrates, 1:12, pos = 4)

gasoline.pcr <- pcr(octane ~ ., data = gasoline,
                    validation = "LOO", ncomp = 1)
RMSEP(gasoline.pcr, estimate = "CV")

gasoline.pcr2 <- pcr(octane ~ ., data = gasoline, ncomp = 1)
X <- gasoline.pcr2$scores
HatM <- X %*% solve(crossprod(X), t(X))
sqrt(mean((gasoline.pcr2$residuals/(1 - diag(HatM)))^2))

sqrt(mean((gasoline.pcr2$residuals/(1 - mean(diag(HatM))))^2))

gasoline.pcr <- pcr(octane ~ ., data = gasoline,
                    validation = "CV", ncomp = 4,
                    segment.type = "interleaved")
RMSEP(gasoline.pcr, estimate = "CV")

gasoline.pls <- plsr(octane ~ ., data = gasoline,
                     validation = "LOO", ncomp = 2,
                     jackknife = TRUE)
n <- length(gasoline$octane)
b.oob <- gasoline.pls$validation$coefficients[,,2,]
bias.est <- (n-1) * (rowMeans(b.oob) - coef(gasoline.pls))
wavelengths <- seq(900, 1700, by = 2)
plot(wavelengths, bias.est, xlab = "wavelength", ylab = "bias",
     type = "h", main = "Jackknife bias estimates")
var.est <- var.jack(gasoline.pls)
lines(wavelengths, var.est, col = "red")

if (boot.present) {
  ## Bootstrap
  B <- 500
  ngas <- nrow(gasoline)
  boot.indices <- 
    matrix(sample(1:ngas, ngas * B, replace = TRUE), ncol = B)
  sort(boot.indices[,1])
  
  npc <- 5
  predictions <- array(NA, c(ngas, npc, B))
  for (i in 1:B) {
    gas.bootpcr <- pcr(octane ~ ., data = gasoline,
                       ncomp = npc, subset = boot.indices[,i])
    oobs <- (1:ngas)[-boot.indices[,i]]
    predictions[oobs,,i] <- 
      predict(gas.bootpcr,
              newdata = gasoline$NIR[oobs,])[,1,]
  }
  diffs <- sweep(predictions, 1, gasoline$octane)
  sqerrors <- apply(diffs^2, c(1,2), mean, na.rm = TRUE)
  sqrt(colMeans(sqerrors))
  
  gas.pcr <- pcr(octane ~ ., data = gasoline, ncomp = npc)
  RMSEP(gas.pcr, intercept = FALSE)
  
  error.632 <- .368 * colMeans(gas.pcr$residuals^2) +
    .632 * colMeans(sqerrors)
  sqrt(error.632)
  
  gas.pcr.cv <- pcr(octane ~ ., data = gasoline, ncomp = npc,
                    validation = "CV")
  gas.pcr.loo <- pcr(octane ~ ., data = gasoline, ncomp = npc,
                     validation = "LOO")
  bp <- barplot(sqrt(error.632),
                ylim = c(0, 1.6), col = "peachpuff")
  lines(bp, sqrt(c(gas.pcr.cv$validation$PRESS) / ngas), 
        col = 2)
  lines(bp, sqrt(c(gas.pcr.loo$validation$PRESS) / ngas), 
        col = 4, lty = 2)
  legend("topright", lty = 1:2, col = c(2, 4), 
         legend = c("CV", "LOO"))
  
  gas.pcr.boot632 <- 
    boot(gasoline, 
         function(x, ind) {
           mod <- pcr(octane ~ ., data = x, subset = ind, ncomp = 4)
           gasoline$octane - predict(mod, newdata = gasoline$NIR, ncomp = 4)
         },
         R = 499)
  dim(gas.pcr.boot632$t)
  
  in.bag <- boot.array(gas.pcr.boot632)
  in.bag[1:10,1]
  
  in.bag <- boot.array(gas.pcr.boot632)
  oob.error <- mean((gas.pcr.boot632$t^2)[in.bag == 0])
  app.error <- MSEP(pcr(octane ~ ., data = gasoline, ncomp = 4),
                    ncomp = 4, intercept = FALSE)
  sqrt(.368 * c(app.error$val) + .632 * oob.error)
  
  B <- 1000
  ngas <- nrow(gasoline)
  boot.indices <- 
    matrix(sample(1:ngas, ngas * B, replace = TRUE), ncol = B)
  npc <- 4
  gas.pcr <- pcr(octane ~ ., data = gasoline, ncomp = npc)
  coefs <- matrix(0, ncol(gasoline$NIR), B)
  for (i in 1:B) {
    gas.bootpcr <- pcr(octane ~ ., data = gasoline,
                       ncomp = npc, subset = boot.indices[,i])
    coefs[,i] <- c(coef(gas.bootpcr))
  }
  matplot(wavelengths, coefs, type = "l",
          lty = 1, col = "gray",
          ylab = "Coefficients", xlab = "Wavelength (nm)")
  abline(h = 0)
  
  coef.stats <- cbind(apply(coefs, 1, quantile, .025),
                      c(coef(gas.pcr)),
                      apply(coefs, 1, quantile, .975))
  matplot(wavelengths, coef.stats, type = "n",
          xlab = "Wavelength (nm)",
          ylab = "Regression coefficient")
  abline(h = 0, col = "gray")
  matlines(wavelengths, coef.stats,
           lty = c(2,1,2), col = c(2,1,2))
  
  gas.pcr.bootCI <-
    boot(gasoline,
         function(x, ind) {
           c(coef(pcr(octane ~ ., data=x, subset = ind)))
         },
         R = 999)
  dim(gas.pcr.bootCI$t)
  
  smallest <- which.min(gas.pcr.bootCI$t0)
  plot(gas.pcr.bootCI, index = smallest)
  
  boot.ci(gas.pcr.bootCI, index = smallest, 
          type = c("perc", "bca"))
}

if (ipred.present) {
  ## Bagging
  odd <- seq(1, length(gasoline$octane), by = 2)
  even <- seq(2, length(gasoline$octane), by = 2)
  gasoline.bagging <- ipredbagg(gasoline$octane[odd],
                                gasoline$NIR[odd,],
                                coob = TRUE)
  gasoline.bagging
  
  gs.baggpreds <- predict(gasoline.bagging, gasoline$NIR[even,])
  resids <- gs.baggpreds - gasoline$octane[even]
  sqrt(mean(resids^2))
}

if (ipred.present) {
  data(prostate, package = "ChemometricsWithRData")
  prost <- prostate[prostate.type != "bph", 1:1000]
  prost.type <- factor(prostate.type[prostate.type != "bph"])
  prost.df <- data.frame(type = prost.type, prost = prost)
  odd <- seq(1, length(prost.type), by = 2)
  even <- seq(2, length(prost.type), by = 2)
  prost.bagging <- bagging(type ~ ., data = prost.df, 
                           subset = odd)
  prost.baggingpred <- predict(prost.bagging, 
                               newdata = prost.df[even,])
  table(prost.type[even], prost.baggingpred)
}

if (randomForest.present) {
  ## Random Forests
  data(wines, package = "ChemometricsWithRData")
  odd <- seq(1, length(vintages), by = 2)
  even <- seq(2, length(vintages), by = 2)
  wines.df <- data.frame(vint = vintages, wines)
  wines.rf <- randomForest(vint ~ ., subset = odd,
                           data = wines.df)
  wines.rf
  
  wines.rf.predict <- predict(wines.rf, 
                              newdata = wines.df[even,])
  sum(wines.rf.predict == wines.df[even,"vint"]) / length(even)
  
  wines.rf <- randomForest(vint ~ ., data = wines.df,
                           importance = TRUE)
  varImpPlot(wines.rf)

  data(gasoline, package = "pls")
  odd <- seq(1, length(gasoline$octane), by = 2)
  even <- seq(2, length(gasoline$octane), by = 2)
  gasoline.rf <- randomForest(gasoline$NIR[odd,], 
                              gasoline$octane[odd],
                              importance = TRUE,
                              xtest = gasoline$NIR[even,],
                              ytest = gasoline$octane[even])
  pl.range <- c(83,90)
  plot(gasoline$octane[odd], gasoline.rf$predicted,
       main = "Training: OOB prediction", xlab = "True",
       ylab = "Predicted", xlim = pl.range, ylim = pl.range)
  abline(0, 1)
  plot(gasoline$octane[even], gasoline.rf$test$predicted,
       main = "Test set prediction", xlab = "True",
       ylab = "Predicted", xlim = pl.range, ylim = pl.range)
  abline(0, 1)
  
  resids <- gasoline.rf$test$predicted - gasoline$octane[even]
  sqrt(mean(resids^2))
  
  rf.imps <- importance(gasoline.rf)
  plot(wavelengths, rf.imps[,1] / max(rf.imps[,1]),
       type = "l", xlab = "Wavelength (nm)",
       ylab = "Importance")
  lines(wavelengths, rf.imps[,2] / max(rf.imps[,2]), col = 2)
  legend("topright", legend = c("Error decrease", "Gini index"),
         col = 1:2, lty = 1)
  
  odd <- seq(1, nrow(prost), by = 2)
  even <- seq(2, nrow(prost), by = 2)
  prost.rf <-
    randomForest(prost[odd,], prost.type[odd],
                 x.test = prost[even,], y.test = prost.type[even])
  prost.rfpred <- predict(prost.rf, newdata = prost[even,])
  table(prost.type[even], prost.rfpred)
}

if (ada.present) {
  ## Boosting
  set.seed(7)
  prost.ada <- ada(type ~ ., data = prost.df, subset = odd)
  prost.adapred <- predict(prost.ada, newdata = prost.df[even,])
  table(prost.type[even], prost.adapred)
  
  prost.ada <- addtest(prost.ada, prost.df[even, -1], prost.type[even])
  plot(prost.ada, test = TRUE)
}
