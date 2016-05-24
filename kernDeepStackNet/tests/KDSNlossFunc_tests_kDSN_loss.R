# Test lossKDSN
library(kernDeepStackNet)
# Check function tune_KDSN (duration is long)
XORdat4Dim <- expand.grid(x1=c(0, 1), x2=c(0, 1), x3=c(0, 1), x4=c(0, 1))
XORdat4Dim <- cbind(y1=as.numeric(xor(xor(xor(XORdat4Dim[, 1], XORdat4Dim[, 2]), XORdat4Dim[, 3]), 
                                      XORdat4Dim[, 4])), XORdat4Dim)
quantEuklid <- quantile(c(dist(robustStandard(XORdat4Dim[, -1]))^2), probs = c(0, 0.5, 1))
# Check one level
as.numeric(lossKDSN (parOpt=c(10, 8, 1), y=matrix(XORdat4Dim[, 1], ncol=1), 
                      X=as.matrix(XORdat4Dim[, -1]), gammaPar=1, seedW=1))
# Check four levels
as.numeric(lossKDSN (parOpt=c(32, 12, 1, 
                    16, 10, 0.1,
                    8, 8, 0.01,
                    4, 6, 0.001), y=matrix(XORdat4Dim[, 1], ncol=1), 
           X=as.matrix(XORdat4Dim[, -1]), gammaPar=1, seedW=1:4))

# Check one level given model
fittedKDSN <- fitKDSN (y=XORdat4Dim[, 1], X=as.matrix(XORdat4Dim[, -1]), 
                       levels=1, Dim=10, 
                       sigma=8, lambda=1, 
                       alpha=0, 
                       info=FALSE, seedW=1, standX=TRUE)
as.numeric(lossKDSNgivenModel (KDSNfit=fittedKDSN, y=matrix(XORdat4Dim[, 1], ncol=1), 
                     X=as.matrix(XORdat4Dim[, -1]), gammaPar=1))

# Check four levels given model
fittedKDSN <- fitKDSN (y=XORdat4Dim[, 1], X=as.matrix(XORdat4Dim[, -1]), 
                       levels=4, Dim=c(32, 16, 8, 4), 
                       sigma=c(12, 10, 8, 6), lambda=c(1, 0.1, 0.01, 0.001), 
                       alpha=rep(0, 4), 
                       info=FALSE, seedW=1:4, standX=TRUE)
as.numeric(lossKDSNgivenModel (KDSNfit=fittedKDSN, y=matrix(XORdat4Dim[, 1], ncol=1), 
                               X=as.matrix(XORdat4Dim[, -1]), gammaPar=1))
