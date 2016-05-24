###################
# Check tuneMboKDSN

library(kernDeepStackNet)
# Check function tune_KDSN (duration is long)
XORdat4Dim <- expand.grid(x1=c(0, 1), x2=c(0, 1), x3=c(0, 1), x4=c(0, 1))
XORdat4Dim <- cbind(y1=as.numeric(xor(xor(xor(XORdat4Dim[, 1], XORdat4Dim[, 2]), XORdat4Dim[, 3]), XORdat4Dim[, 4])), 
                    XORdat4Dim)
tunedKDSN <- tuneMboKDSN (y=XORdat4Dim[, 1], X=as.matrix(XORdat4Dim[, -1]),  
                        maxLevels=2, gammaPar=1, nStepMult=1, designMult=5, fineTuneIt=10, dimMax=round(dim(as.matrix(XORdat4Dim[, -1]))[1]/4))
library(pROC)
cat(all.equal(as.numeric(auc(XORdat4Dim[, 1], c(predict(tunedKDSN, newx=as.matrix(XORdat4Dim[, -1]))))), 1))

# Check if gcv score is available
stopifnot(!is.null(attr(tunedKDSN, which="GCV")))
