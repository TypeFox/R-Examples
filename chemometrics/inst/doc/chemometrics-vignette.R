### R code from vignette source 'chemometrics-vignette.rnw'

###################################################
### code chunk number 1: chemometrics-vignette.rnw:217-220
###################################################
  ## set the prompt to "> " and the continuation to "  "
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 2: chemometrics-vignette.rnw:258-259
###################################################
  library(chemometrics)


###################################################
### code chunk number 3: chemometrics-vignette.rnw:332-334
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 4: chemometrics-vignette.rnw:428-430
###################################################
  data(wine,package="gclus")
  X <- scale(wine[,2:14])       # first column: wine type


###################################################
### code chunk number 5: chemometrics-vignette.rnw:440-441 (eval = FALSE)
###################################################
##   res <- pcaCV(X)


###################################################
### code chunk number 6: 01-pcaCV
###################################################
  par(mfrow=c(1,2))
  par(mar=c(5, 4, 1, 1))
  pcaCV(X)
  par(mar=c(7.5, 4, 1, 1))
  pcaVarexpl(X, a=5, las=2)


###################################################
### code chunk number 7: chemometrics-vignette.rnw:477-478
###################################################
  X_nipals <- nipals(X, a=5)


###################################################
### code chunk number 8: chemometrics-vignette.rnw:535-536
###################################################
  X_nipals <- nipals(X, a=5, it=160)


###################################################
### code chunk number 9: chemometrics-vignette.rnw:559-560 (eval = FALSE)
###################################################
##   res <- pcaVarexpl(X, a=5)


###################################################
### code chunk number 10: chemometrics-vignette.rnw:613-615
###################################################
  X_nipals <- list(scores=X_nipals$T, loadings=X_nipals$P, 
    sdev=apply(X_nipals$T, 2, sd))


###################################################
### code chunk number 11: chemometrics-vignette.rnw:617-618 (eval = FALSE)
###################################################
##   res <- pcaDiagplot(X, X.pca=X_nipals, a=5)


###################################################
### code chunk number 12: 03-pcaDiagplot-nipals
###################################################
  par(mar=c(5,4,1,1))
  res <- pcaDiagplot(X, X.pca=X_nipals, a=5)


###################################################
### code chunk number 13: chemometrics-vignette.rnw:636-641 (eval = FALSE)
###################################################
##   X.rob <- scale(wine[,2:14], center = apply(wine[,2:14], 2, median), 
##     scale = apply(wine[,2:14], 2, mad))
##   library(pcaPP)                 # robust PCA based on projection pursuit
##   X.grid <- PCAgrid(X.rob,k=5,scale=mad)
##   res1 <- pcaDiagplot(X.rob,X.grid,a=5)


###################################################
### code chunk number 14: 04-pcaDiagplot-robust
###################################################
  X.rob <- scale(wine[,2:14], center = apply(wine[,2:14], 2, median), 
    scale = apply(wine[,2:14], 2, mad))
  library(pcaPP)                 # robust PCA based on projection pursuit
  X.grid <- PCAgrid(X.rob,k=5,scale=mad)
  par(mar=c(5,4,1,1))
  res1 <- pcaDiagplot(X.rob,X.grid,a=5)


###################################################
### code chunk number 15: chemometrics-vignette.rnw:663-680
###################################################
# diagnostics with robust PCA
  orth1 <- dimnames(X)[[1]][(res1$ODist > res1$critOD) * (res1$SDist < res1$critSD) == 1] # orthogonal outliers
  good1 <- dimnames(X)[[1]][(res1$ODist < res1$critOD) * (res1$SDist > res1$critSD) == 1] # good leverage points
  bad1 <- dimnames(X)[[1]][(res1$ODist > res1$critOD) * (res1$SDist > res1$critSD) == 1] # bad leverage points
# diagnostics with NIPALS
  orth2 <- dimnames(X)[[1]][(res$ODist > res$critOD) * (res$SDist < res$critSD) == 1] # orthogonal outliers
  good2 <- dimnames(X)[[1]][(res$ODist < res$critOD) * (res$SDist > res$critSD) == 1] # good leverage points
  bad2 <- dimnames(X)[[1]][(res$ODist > res$critOD) * (res$SDist > res$critSD) == 1] # bad leverage points
# print:
  cat("DIAGNOSTICS WITH NIPALS:", "\n",
      "orthogonal outliers: ", orth2, "\n",
      "good leverage points:", good2, "\n",
      "bad leverage points: ", bad2, "\n", sep=" ")
  cat("DIAGNOSTICS WITH ROBUST PCA:", "\n",
      "orthogonal outliers: ", orth1, "\n",
      "good leverage points:", good1, "\n",
      "bad leverage points: ", bad1, "\n", sep=" ")


###################################################
### code chunk number 16: chemometrics-vignette.rnw:691-693
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 17: chemometrics-vignette.rnw:742-747
###################################################
  data(NIR)
  X <- NIR$xNIR
  Y <- NIR$yGlcEtOH
  namesX <- names(X)
  attach(X)


###################################################
### code chunk number 18: chemometrics-vignette.rnw:749-754 (eval = FALSE)
###################################################
##   data(NIR)
##   X <- NIR$xNIR
##   Y <- NIR$yGlcEtOH
##   namesX <- names(X)
##   attach(X)


###################################################
### code chunk number 19: chemometrics-vignette.rnw:821-823
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 20: chemometrics-vignette.rnw:864-866
###################################################
  y <- Y[,1]
  NIR.Glc <- data.frame(X=X, y=y)


###################################################
### code chunk number 21: chemometrics-vignette.rnw:868-869 (eval = FALSE)
###################################################
##   res.step <- stepwise(y~., data=NIR.Glc)


###################################################
### code chunk number 22: chemometrics-vignette.rnw:878-882 (eval = FALSE)
###################################################
##   steps    <- length(res.step$usedTime)
##   seconds  <- res.step$usedTime[steps]
##   varnbr   <- sum(finalvars <- res.step$models[steps,])
##   varnames <- namesX[as.logical(finalvars)]


###################################################
### code chunk number 23: chemometrics-vignette.rnw:885-890
###################################################
  steps <- 22
  seconds <- 15
  varnbr <- 16
  varnames <- c('X1115.0','X1185.0','X1215.0','X1385.0','X1420.0','X1500.0','X1565.0','X1585.0',
                'X1690.0','X1715.0','X1720.0','X1815.0','X1995.0','X2070.0','X2100.0','X2195.0')


###################################################
### code chunk number 24: chemometrics-vignette.rnw:893-896
###################################################
  cat(" steps needed:        ", steps, "\n",
      "seconds needed:      ", seconds, "\n",
      "number of variables: ", varnbr, "\n", sep=" ")


###################################################
### code chunk number 25: 01-BIC (eval = FALSE)
###################################################
##   # produce Figure 5
##   modelsize <- apply(res.step$models, 1, sum)
##   plot(modelsize, res.step$bic, type="b", pch=20,
##     main="BIC during stepwise regression",
##     xlab="model size", ylab="BIC value")


###################################################
### code chunk number 26: 01-BIC
###################################################
  # produce Figure 4.1
  modelsize <- c(1,2,3,4,5,6,7,8,9,10,11,10,11,12,11,12,13,14,15,16,15,16)
  res.step.bic <- c(1292.227,1231.984,1216.206,1182.616,1176.313,1172.719,1136.370,1075.659,
    1064.525,1060.948,1057.010,1053.532,1050.034,1047.147,1042.120,1039.618,
    1032.290,1023.446,1022.557,1022.301,1018.909,1017.506)
  plot(modelsize, res.step.bic, type="b", pch=20,
    main="BIC during stepwise regression",
    xlab="model size", ylab="BIC value")


###################################################
### code chunk number 27: chemometrics-vignette.rnw:937-940
###################################################
  finalmodel <- as.formula(paste("y~", paste(varnames, collapse = "+"),
    sep = ""))
  res.stepcv <- lmCV(finalmodel, data=NIR.Glc, repl=100, segments=4)


###################################################
### code chunk number 28: 02-lmCV
###################################################
  par(mfrow=c(1,2))
  plot(rep(y,100), res.stepcv$predicted, pch=".",
    main="Predictions from repeated CV",
    xlab="Measured y", ylab="Predicted y")
  abline(coef=c(0,1))
  points(y, apply(res.stepcv$predicted, 1, mean), pch=20)
  plot(res.stepcv$SEP, main="Standard Error of Prediction",
    xlab="Number of repetitions", ylab="SEP")
  abline(h=res.stepcv$SEPm)
  abline(h=median(res.stepcv$SEP), lty=2)


###################################################
### code chunk number 29: chemometrics-vignette.rnw:976-979
###################################################
    cat(" SEP mean:               ", round(res.stepcv$SEPm, 3), "\n",
      "SEP median:             ", round(median(res.stepcv$SEP), 3), "\n",
      "SEP standard deviation: ", round(sd(res.stepcv$SEP), 3), "\n", sep=" ")


###################################################
### code chunk number 30: chemometrics-vignette.rnw:997-999
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 31: chemometrics-vignette.rnw:1057-1059 (eval = FALSE)
###################################################
##   res.pcr <- mvr_dcv(y~., data=NIR.Glc, ncomp=40, method = "svdpc",
##     repl = 100)


###################################################
### code chunk number 32: 03-comppcr (eval = FALSE)
###################################################
##   plotcompmvr(res.pcr)


###################################################
### code chunk number 33: 04-seppcr (eval = FALSE)
###################################################
##   optpcr <- res.pcr$afinal
##   plotSEPmvr(res.pcr, optcomp=optpcr, y=y, X=X, method="svdpc")


###################################################
### code chunk number 34: chemometrics-vignette.rnw:1132-1133
###################################################
  optpcr <- 31


###################################################
### code chunk number 35: 05-predpcr (eval = FALSE)
###################################################
##   plotpredmvr(res.pcr, optcomp=optpcr, y=y, X=X, method="svdpc")


###################################################
### code chunk number 36: 06-respcr (eval = FALSE)
###################################################
##   plotresmvr(res.pcr, optcomp=optpcr, y=y, X=X, method="svdpc")


###################################################
### code chunk number 37: chemometrics-vignette.rnw:1182-1187
###################################################
    #SEP <- apply(res.pcr$resopt[,1,], 2, sd)
    cat(" optimal number of PCs:  ", optpcr, "\n",
         "SEP mean:               7.433 \n",                      #round(mean(SEP), 4)
         "SEP median:             7.449 \n",                     #round(median(SEP), 4)
         "SEP standard deviation: 0.373 \n", sep=" ")    # round(sd(SEP), 4)


###################################################
### code chunk number 38: chemometrics-vignette.rnw:1205-1207
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 39: chemometrics-vignette.rnw:1462-1464 (eval = FALSE)
###################################################
##   res.pls <- mvr_dcv(y~., data=NIR.Glc, ncomp=40, method = "simpls",
##     repl = 100)


###################################################
### code chunk number 40: 07-comppls (eval = FALSE)
###################################################
##   plotcompmvr(res.pls)


###################################################
### code chunk number 41: 08-seppls (eval = FALSE)
###################################################
##   optpls <- res.pls$afinal
##   plotSEPmvr(res.pls, optcomp=optpls, y=y, X=X, method="simpls")


###################################################
### code chunk number 42: chemometrics-vignette.rnw:1479-1480
###################################################
  optpls <- 14


###################################################
### code chunk number 43: 09-predpls (eval = FALSE)
###################################################
##   plotpredmvr(res.pls, optcomp=optpls, y=y, X=X , method="simpls")


###################################################
### code chunk number 44: 10-respls (eval = FALSE)
###################################################
##   plotresmvr(res.pls, optcomp=optpls, y=y, X=X, method="simpls")


###################################################
### code chunk number 45: chemometrics-vignette.rnw:1542-1547
###################################################
    #SEP <- apply(res.pls$resopt[,1,], 2, sd)
    cat(" optimal number of PCs:  ", optpls, "\n",
         "SEP mean:                6.489 \n", #round(mean(SEP), 4)
         "SEP median:              6.491 \n", #round(median(SEP), 4)
         "SEP standard deviation:  0.408 \n", sep=" ") #round(sd(SEP), 4)


###################################################
### code chunk number 46: chemometrics-vignette.rnw:1561-1562 (eval = FALSE)
###################################################
##   res.pls1nipals <- pls1_nipals(X, y, a = res.pls$afinal)


###################################################
### code chunk number 47: chemometrics-vignette.rnw:1569-1571 (eval = FALSE)
###################################################
##   Ysc <- scale(Y)
##   res.pls2nipals <- pls2_nipals(X, Ysc, a = 9)


###################################################
### code chunk number 48: chemometrics-vignette.rnw:1576-1578 (eval = FALSE)
###################################################
##   a <- min(dim(X)[1], dim(X)[2], dim(Y)[2])
##   res.plseigen <- pls_eigen(X, Ysc, a = a)


###################################################
### code chunk number 49: chemometrics-vignette.rnw:1594-1596
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 50: 11-prmcv (eval = FALSE)
###################################################
##   res.prmcv <- prm_cv(X, y, a = 40, opt = "median")


###################################################
### code chunk number 51: chemometrics-vignette.rnw:1775-1781
###################################################
  oc <- 20
  cat(" optimal number of PCs: 20 \n",  #oc <- res.prmcv$optcomp
       "                         classic   20% trimmed \n",
       "SEP mean:                5.454      5.094 \n",          #round(mean(res.prmcv$SEPj[,oc]), 4), round(mean(res.prmcv$SEPtrimj[,oc]), 4)
       "SEP median:              5.488      5.070 \n",          # round(median(res.prmcv$SEPj[,oc]), 4), round(median(res.prmcv$SEPtrimj[,oc]), 4)
       "SEP standard deviation:  0.627      1.289 \n", sep=" ") #round(sd(res.prmcv$SEPj[,oc]), 4), round(sd(res.prmcv$SEPtrimj[,oc]), 4)


###################################################
### code chunk number 52: 12-plotprm (eval = FALSE)
###################################################
##   plotprm(res.prmcv, y)


###################################################
### code chunk number 53: chemometrics-vignette.rnw:1807-1808 (eval = FALSE)
###################################################
##   prm(X, y, a = res.prmcv$optcomp, opt = "l1m", usesvd = TRUE)


###################################################
### code chunk number 54: 11a-prmdcv (eval = FALSE)
###################################################
##   res.prmdcv <- prm_dcv(X, y, a = 40, opt = "median",repl=20)


###################################################
### code chunk number 55: 07a-compprmdcv (eval = FALSE)
###################################################
##   plotcompprm(res.prmdcv)


###################################################
### code chunk number 56: 08a-sepprmdcv (eval = FALSE)
###################################################
##   plotSEPprm(res.prmdcv,res.prmdcv$afinal,y,X)


###################################################
### code chunk number 57: 08b-predprmdcv (eval = FALSE)
###################################################
##   plotpredprm(res.prmdcv,res.prmdcv$afinal,y,X)


###################################################
### code chunk number 58: 08c-resprmdcv (eval = FALSE)
###################################################
## plotresprm(res.prmdcv,res.prmdcv$afinal,y,X)


###################################################
### code chunk number 59: chemometrics-vignette.rnw:1908-1910
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 60: 13-plotridge
###################################################
  res.plotridge <- plotRidge(y~., data=NIR.Glc, lambda=seq(0.5,10,by=0.05))


###################################################
### code chunk number 61: 14-ridgecv (eval = FALSE)
###################################################
##   res.ridge <- ridgeCV(y~., data=NIR.Glc, lambdaopt=res.plotridge$lambdaopt,
##     repl=100)


###################################################
### code chunk number 62: chemometrics-vignette.rnw:2055-2057
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 63: 15-lassocv (eval = FALSE)
###################################################
##   res.lasso <- lassoCV(y~., data = NIR.Glc, fraction = seq(0, 1, by = 0.05), 
##     legpos="top")


###################################################
### code chunk number 64: 16-lassocoef (eval = FALSE)
###################################################
##   res.lassocoef <- lassocoef(y~., data=NIR.Glc, sopt=res.lasso$sopt)


###################################################
### code chunk number 65: chemometrics-vignette.rnw:2175-2177
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 66: chemometrics-vignette.rnw:2200-2225
###################################################
  SEP.stepwise <- 4.61 #res.stepcv$SEPm
  SEP20.stepwise <- 4.40 #mean(apply(resid.trim.stepwise, 2, sd))
  SEP.pcr <- 7.43 #mean(apply(res.pcr$resopt[,1,], 2, sd))
  SEP20.pcr <- 7.37 #mean(apply(resid.trim.pcr, 2, sd))
  SEP.pls <- 6.49 #mean(apply(res.pls$resopt[,1,], 2, sd))
  SEP20.pls <- 6.52 #mean(apply(resid.trim.pls, 2, sd))
  SEP.prm <- 5.40 #res.prmcv$SEPall[res.prmcv$optcomp]
  SEP20.prm <- 4.95 #res.prmcv$SEPop
  SEP.prmdcv <- 5.95 #res.prmdcv$SEPall[res.prmcv$optcomp]
  SEP20.prmdcv <- 5.86 #res.prmcdv$SEPop
  SEP.ridge <- 5.37 #res.ridge$SEPm
  SEP20.ridge <- 5.32 #mean(apply(resid.trim.ridge, 2, sd))
  SEP.lasso <- 6.48 #res.lasso$SEP[(res.lasso$fraction==res.lasso$sopt)]
  SEP20.lasso <- 5.89 #compute trimmed SEP within lassoCV function (not implemented)
  cat(" PREDICTION PERFORMANCE",
    "\n Method      SEP     SEP20%     Nr. of Variables / Components ",
      "\n                                    ",
    "\n stepwise    ", round(SEP.stepwise, 2), "    ", round(SEP20.stepwise, 2), "        ", varnbr, " variables",
    "\n PCR         ", round(SEP.pcr, 2), "    ", round(SEP20.pcr, 2), "       ", optpcr, " components",
    "\n PLS         ", round(SEP.pls, 2), "    ", round(SEP20.pls, 2), "       ", optpls, " components",
    "\n PRM-CV      ", round(SEP.prm, 2), "     ", round(SEP20.prm, 2), "       ", oc, " components",
    "\n PRM-DCV     ", round(SEP.prmdcv, 2), "    ", round(SEP20.prmdcv, 2), "       ", oc, " components",
    "\n Ridge       ", round(SEP.ridge, 2), "    ", round(SEP20.ridge, 2), "      ", dim(X)[2], " variables",
    "\n Lasso       ", round(SEP.lasso, 2), "    ", round(SEP20.lasso, 2), "      110 variables \n", #res.lassocoef$numb.nonzero
    sep="")


###################################################
### code chunk number 67: chemometrics-vignette.rnw:2232-2241
###################################################
  cat(" COMPUTATION TIMES \n",
    "Method      Algorithm     Time needed \n \n",
    "stepwise    100 x RCV       0min 44sec \n",
    "PCR         100 x RDCV      2min 56sec \n",
    "PLS         100 x RDCV      2min 27sec \n",
    "PRM         single CV       4min 15sec \n",
    "PRM         20 x RDCV     241min 15sec \n",
    "Ridge       100 x RCV       1min 40sec \n",
    "Lasso       single CV       0min 33sec \n")


###################################################
### code chunk number 68: chemometrics-vignette.rnw:2258-2260
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 69: chemometrics-vignette.rnw:2332-2338
###################################################
  library(gclus)
  data(wine)
  X <- data.frame(scale(wine[,2:14]))   # data without class information
  grp <- as.factor(wine[,1])            # class information
  wine <- data.frame(X=X, grp=grp)
  train <- sample(1:length(grp), round(2/3*length(grp)))


###################################################
### code chunk number 70: chemometrics-vignette.rnw:2348-2350
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 71: 01-knneval
###################################################
  knneval <- knnEval(X, grp, train, knnvec=seq(1,30), legpos="topright")


###################################################
### code chunk number 72: chemometrics-vignette.rnw:2393-2398
###################################################
  indmin <- which.min(knneval$cvMean)
  thresh <- knneval$cvMean[indmin] + knneval$cvSe[indmin]
  fvec <- (knneval$cvMean < thresh)
  indopt <- min((1:indmin)[fvec[1:indmin]])
  opt <- knneval$knnvec[indopt]


###################################################
### code chunk number 73: chemometrics-vignette.rnw:2408-2411
###################################################
  cat(" optimal number of nearest neighbors:", "k =", opt, "\n",
      "test error at optimum: ", round(knneval$testerr[indopt],4), "\n",
      "CV error threshold:    ", round(thresh,4), "\n", sep=" ")


###################################################
### code chunk number 74: 02-repknn (eval = FALSE)
###################################################
##   res <- array(dim=c(100,6))
##   colnames(res) <- c("k","trainerr","testerr","cvmean","cvse","threshold")
##   tt <- proc.time()
##   for (i in 1:100) {
##     train <- sample(1:length(grp), round(2/3*length(grp)))
##     knneval <- knnEval(X, grp, train, knnvec=seq(1,30), plotit=FALSE)
##     indmin <- which.min(knneval$cvMean)
##     res[i,6] <- knneval$cvMean[indmin] + knneval$cvSe[indmin]
##     fvec <- (knneval$cvMean < res[i,6])
##     indopt <- min((1:indmin)[fvec[1:indmin]])
##     res[i,1] <- knneval$knnvec[indopt]
##     res[i,2] <- knneval$trainerr[indopt]
##     res[i,3] <- knneval$testerr[indopt]
##     res[i,4] <- knneval$cvMean[indopt]
##     res[i,5] <- knneval$cvSe[indopt]
##   }
##   tt <- proc.time() - tt
##   plot(table(res[,1]), xlab="Number of Nearest Neighbors", ylab="Frequency")


###################################################
### code chunk number 75: chemometrics-vignette.rnw:2456-2464
###################################################
  mins <- 263.52/60; leftsec <- (mins - floor(mins))*60
  cat("                   median    sd", "\n",
       " training error   0.0169   0.0138 \n", #round(median(res[,2], na.rm=TRUE),4), round(sd(res[,2], na.rm=TRUE),4), "\n",
       " test error       0.05     0.0228 \n", #round(median(res[,3], na.rm=TRUE),4), " ", round(sd(res[,3], na.rm=TRUE),4), "\n",
       " CV error mean    0.0258   0.0136 \n", #round(median(res[,4], na.rm=TRUE),4), round(sd(res[,4], na.rm=TRUE),4), "\n",
       " (computing time for 100 repetitions: 4 min 24 sec) \n", sep=" ") #, floor(mins), "min", round(leftsec,0), "sec) \n", sep=" ")
  #resknn <- res
  timeknn <- round(mins,0)


###################################################
### code chunk number 76: chemometrics-vignette.rnw:2471-2472 (eval = FALSE)
###################################################
##   pred <- knn(X[train,], X[-train,], cl=grp[train], k = 1)


###################################################
### code chunk number 77: chemometrics-vignette.rnw:2479-2481
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 78: chemometrics-vignette.rnw:2545-2549 (eval = FALSE)
###################################################
##   # library(rpart)
##   tree <- rpart::rpart(grp~., data=wine, method="class")
##   plot(tree, main="Full Tree")
##   text(tree)


###################################################
### code chunk number 79: chemometrics-vignette.rnw:2560-2561 (eval = FALSE)
###################################################
##   treeeval <- treeEval(X, grp, train, cp=seq(0.45,0.05,by=-0.05))


###################################################
### code chunk number 80: 03-treeeval
###################################################
  tree <- rpart(grp~., data=wine, method="class")
  par(mfrow=c(1,2))
  par(mar=c(1,0,4,2))
  plot(tree, main="Full Tree")
  text(tree)
  par(mar=c(5,4,1,2))
  treeeval <- treeEval(X, grp, train, cp=seq(0.45,0.05,by=-0.05), legpos="topright")


###################################################
### code chunk number 81: chemometrics-vignette.rnw:2598-2606
###################################################
  indmin <- which.min(treeeval$cvMean)
  thresh <- treeeval$cvMean[indmin] + treeeval$cvSe[indmin]
  fvec <- (treeeval$cvMean < thresh)
  indopt <- min((1:indmin)[fvec[1:indmin]])
  opt <- treeeval$cp[indopt]
  cat(" optimal tree complexity:", "cp =", opt, "\n",
      "test error at optimum: ", round(treeeval$testerr[indopt],4), "\n",
      "CV error threshold:    ", round(thresh,4), "\n", sep=" ")


###################################################
### code chunk number 82: chemometrics-vignette.rnw:2613-2631 (eval = FALSE)
###################################################
##   res <- array(dim=c(100,6))
##   colnames(res) <- c("cp","trainerr","testerr","cvmean","cvse","threshold")
##   tt <- proc.time()
##   for (i in 1:100) {
##     print(i)
##     train <- sample(1:length(grp), round(2/3*length(grp)))
##     treeeval <- treeEval(X, grp, train, cp=seq(0.45,0.05,by=-0.05), plotit=FALSE)
##     indmin <- which.min(treeeval$cvMean)
##     res[i,6] <- treeeval$cvMean[indmin] + treeeval$cvSe[indmin]
##     fvec <- (treeeval$cvMean < res[i,6])
##     indopt <- min((1:indmin)[fvec[1:indmin]])
##     res[i,1] <- treeeval$cp[indopt]
##     res[i,2] <- treeeval$trainerr[indopt]
##     res[i,3] <- treeeval$testerr[indopt]
##     res[i,4] <- treeeval$cvMean[indopt]
##     res[i,5] <- treeeval$cvSe[indopt]
##   }
##   tt <- proc.time() - tt


###################################################
### code chunk number 83: chemometrics-vignette.rnw:2637-2645
###################################################
  mins <- 450.81/60; leftsec <- (mins - floor(mins))*60
  cat("                   median     sd \n",
       " training error   0.0763   0.0325 \n", #round(median(res[,2], na.rm=TRUE),4), round(sd(res[,2], na.rm=TRUE),4), "\n",
       " test error       0.15     0.045  \n", #round(median(res[,3], na.rm=TRUE),4), " ", round(sd(res[,3], na.rm=TRUE),4), "\n",
       " CV error mean    0.1432   0.04   \n", #round(median(res[,4], na.rm=TRUE),4), round(sd(res[,4], na.rm=TRUE),4), "\n",
       " (computing time for 100 repetitions: 7 min 31 sec) \n", sep=" ") #floor(mins), "min", round(leftsec,0), "sec) \n", sep=" ")
  #restree <- res
  timetree <- round(mins,0)


###################################################
### code chunk number 84: chemometrics-vignette.rnw:2654-2657 (eval = FALSE)
###################################################
##   opttree <- prune(tree, cp=0.3)
##   plot(opttree, main="Optimal Tree")
##   text(opttree)


###################################################
### code chunk number 85: 04-prunetree (eval = FALSE)
###################################################
##   par(mfrow=c(1,2))
##   par(mar=c(5,4,1,0))
##   plot(table(res[,1]), xlab="Complexity Parameters", ylab="Frequency")
##   opttree <- prune(tree, cp=0.3)
##   par(mar=c(2,4,4,0))
##   plot(opttree, main="Optimal Tree")
##   text(opttree)


###################################################
### code chunk number 86: chemometrics-vignette.rnw:2675-2677
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 87: chemometrics-vignette.rnw:2743-2746 (eval = FALSE)
###################################################
##   weightsel <- c(0,0.01,0.1,0.15,0.2,0.3,0.5,1)
##   nneteval <- nnetEval(X, grp, train, decay=weightsel, size=5)
##   nneteval <- nnetEval(X, grp, train, decay=0.2, size=seq(5,30,by=5))


###################################################
### code chunk number 88: chemometrics-vignette.rnw:2759-2778 (eval = FALSE)
###################################################
##   res <- array(dim=c(100,6))
##   colnames(res) <-
##     c("weight","trainerr","testerr","cvmean","cvse","threshold")
##   tt <- proc.time()
##   for (i in 1:100) {
##     print(i)
##     train <- sample(1:length(grp), round(2/3*length(grp)))
##     nneteval <- nnetEval(X, grp, train, decay=weightsel, size=5, plotit=FALSE)
##     indmin <- which.min(nneteval$cvMean)
##     res[i,6] <- nneteval$cvMean[indmin] + nneteval$cvSe[indmin]
##     fvec <- (nneteval$cvMean < res[i,6])
##     indopt <- min((1:indmin)[fvec[1:indmin]])
##     res[i,1] <- nneteval$decay[indopt]
##     res[i,2] <- nneteval$trainerr[indopt]
##     res[i,3] <- nneteval$testerr[indopt]
##     res[i,4] <- nneteval$cvMean[indopt]
##     res[i,5] <- nneteval$cvSe[indopt]
##   }
##   tt <- proc.time() - tt


###################################################
### code chunk number 89: chemometrics-vignette.rnw:2784-2792
###################################################
  mins <- 614.66/60; leftsec <- (mins - floor(mins))*60
  cat("                   median     sd", "\n",
       " training error   0        0.0013 \n", #round(median(res[,2], na.rm=TRUE),4), "    ", round(sd(res[,2], na.rm=TRUE),4), "\n",
       " test error       0.0169   0.0167 \n", #round(median(res[,3], na.rm=TRUE),4), round(sd(res[,3], na.rm=TRUE),4), "\n",
       " CV error mean    0.0167   0.0091 \n", #round(median(res[,4], na.rm=TRUE),4), round(sd(res[,4], na.rm=TRUE),4), "\n",
       " (computing time for 100 repetitions: 10 min 15 sec) \n", sep=" ") #floor(mins), "min", round(leftsec,0), "sec) \n", sep=" ")
  #resnnet <- res
  timennet <- round(mins,0)


###################################################
### code chunk number 90: 05-nnetfreq (eval = FALSE)
###################################################
##   plot(table(res[,1]), xlab="weights", ylab="Frequency")


###################################################
### code chunk number 91: 06-nneteval (eval = FALSE)
###################################################
##   par(mfrow=c(1,2))
##   nneteval <- nnetEval(wine, grp, train, size=5, legpos="topright",
##     decay =  c(0,0.01,0.1,0.15,0.2,0.3,0.5,1))
##   nneteval <- nnetEval(wine, grp, train, decay=0.01, legpos="topright",
##     size = c(5,10,15,20,25,30,40))


###################################################
### code chunk number 92: chemometrics-vignette.rnw:2830-2834 (eval = FALSE)
###################################################
##   rule <- nnet(X[train,], class.ind(grp[train]), size=5, entropy=TRUE,
##     decay=0.01)
##   pred <- predict(rule, X[-train,])  # predicted probabilities for test set
##   pred <- apply(pred,1,which.max)    # predicted groups for test set


###################################################
### code chunk number 93: chemometrics-vignette.rnw:2841-2843
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 94: chemometrics-vignette.rnw:2954-2957 (eval = FALSE)
###################################################
##   gamsel <- c(0.001,0.01,0.02,0.05,0.1,0.15,0.2,0.5,1)
##   svmeval <- svmEval(X, grp, train, gamvec=gamsel, kernel="radial",
##     legpos="top")


###################################################
### code chunk number 95: chemometrics-vignette.rnw:2969-2988 (eval = FALSE)
###################################################
##   repl <- 100
##   gamsel <- c(0.001,0.01,0.02,0.05,0.1,0.15,0.2,0.5,1)
##   res <- array(dim=c(repl,6))
##   colnames(res) <- c("gamma","trainerr","testerr","cvmean","cvse","threshold")
##   tt <- proc.time()
##   for (i in 1:repl) {
##     train <- sample(1:length(grp), round(2/3*length(grp)))
##     se <- svmEval(X, grp, train, gamvec=gamsel, kernel="radial", plotit=FALSE)
##     indmin <- which.min(se$cvMean)
##     res[i,6] <- se$cvMean[indmin] + se$cvSe[indmin]
##     fvec <- (se$cvMean < res[i,6])
##     indopt <- min((1:indmin)[fvec[1:indmin]])
##     res[i,1] <- se$gamvec[indopt]
##     res[i,2] <- se$trainerr[indopt]
##     res[i,3] <- se$testerr[indopt]
##     res[i,4] <- se$cvMean[indopt]
##     res[i,5] <- se$cvSe[indopt]
##   }
##   tt <- proc.time() - tt


###################################################
### code chunk number 96: 08-svmeval (eval = FALSE)
###################################################
##   par(mfrow=c(1,2))
##   gamsel <- c(0.001,0.01,0.02,0.05,0.1,0.15,0.2,0.5,1)
##   svmeval <- svmEval(X, grp, train, gamvec=gamsel, kernel="radial",
##     legpos="top")
##   plot(table(res[,1]), xlab="Gamma", ylab="Frequency")


###################################################
### code chunk number 97: chemometrics-vignette.rnw:3006-3014
###################################################
  mins <- 435.87/60; leftsec <- (mins - floor(mins))*60
  cat("                   median     sd", "\n",
       " training error   0.0085   0.0061 \n", #round(median(res[,2], na.rm=TRUE),4), round(sd(res[,2], na.rm=TRUE),4), "\n",
       " test error       0.0167   0.0145 \n", #round(median(res[,3], na.rm=TRUE),4), round(sd(res[,3], na.rm=TRUE),4), "\n",
       " CV error mean    0.0167   0.0081 \n", #round(median(res[,4], na.rm=TRUE),4), round(sd(res[,4], na.rm=TRUE),4), "\n",
       " (computing time for 100 repetitions: 7 min 16 sec) \n", sep=" ") #floor(mins), "min", round(leftsec,0), "sec) \n", sep=" ")
 #ressvm <- res
  timesvm <- round(mins,0)


###################################################
### code chunk number 98: chemometrics-vignette.rnw:3022-3024 (eval = FALSE)
###################################################
##   rule <- svm(X[train,],grp[train], kernel="radial", gamma=0.01)
##   pred <- predict(svmres, X[-train,])


###################################################
### code chunk number 99: chemometrics-vignette.rnw:3031-3033
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 100: 09-classcomp (eval = FALSE)
###################################################
##   par(mar=c(2,4,2,1))
##   resdat <- data.frame(kNN=resknn[,3],Tree=restree[,3],ANN=resnnet[,3],SVM=ressvm[,3])
##   boxplot(resdat,ylab="Test error", las=1)


###################################################
### code chunk number 101: chemometrics-vignette.rnw:3063-3068
###################################################
  cat("Method   Median test error   Computing time \n \n",
      "kNN            0.05              4 min \n", #round(median(resknn[,3], na.rm=TRUE),3),  "              ", timeknn,  " min \n",
      "Tree           0.15              8 min \n", #round(median(restree[,3], na.rm=TRUE),3), "              ", timetree, " min \n",
      "ANN            0.017            10 min \n", #round(median(resnnet[,3], na.rm=TRUE),3), "            ", timennet,   " min \n",
      "SVM            0.017             7 min \n", sep="") #round(median(ressvm[,3], na.rm=TRUE),3),  "             ", timesvm,   " min \n", sep="")


###################################################
### code chunk number 102: chemometrics-vignette.rnw:3192-3194
###################################################
  options(prompt="> ", continue="  ")
  options(width=70)


###################################################
### code chunk number 103: chemometrics-vignette.rnw:3314-3315 (eval = FALSE)
###################################################
##   X_alr <- alr(X, divisorvar=2)


###################################################
### code chunk number 104: chemometrics-vignette.rnw:3337-3338 (eval = FALSE)
###################################################
##   X_clr <- clr(X)


###################################################
### code chunk number 105: chemometrics-vignette.rnw:3371-3372 (eval = FALSE)
###################################################
##   X_ilr <- ilr(X)


