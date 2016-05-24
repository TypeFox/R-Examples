library(caret)
library(kernDeepStackNet)

# Construct test cases with linear model
set.seed(0)
treeTestInd <- createDataPartition(y=trees$Volume, p = 0.8, list = TRUE, times=10)

# Fit linear model
err <- vector("numeric", 10)
for(i in 1:10) {
  lmPart <- lm(Volume~Height+Girth, data=trees[treeTestInd[[i]], ])
  preds <- predict(lmPart, newdata=trees[-treeTestInd[[i]], ])
  err[i] <- sqrt(mean((trees[-treeTestInd[[i]], "Volume"]-preds)^2))
}
mean(err)

# Fit with KDSN and three levels
tempMat <- robustStandard(as.matrix(trees[treeTestInd[[1]], ]))
tempMat <- dist(tempMat)
tempVec <- c(tempMat^2)
quantEuklid <- quantile(tempVec, probs = c(0.25, 0.75))
errKDSN <- vector("numeric", 1)
Level <- 3
for(i in 1:10) {
  KDSNpart <- fitKDSN(y=trees[treeTestInd[[i]], "Volume"], X=as.matrix(trees[treeTestInd[[i]], -3]), 
                      levels=Level, Dim=round(seq(dim(trees)[1], sqrt(dim(trees)[1]), length.out=Level)), 
                      sigma=seq(quantEuklid[1], quantEuklid[2], length.out=Level), lambda=seq(10^-1, 10^-10, length.out=Level), 
                      alpha=rep(0, Level), 
                      info=FALSE, seedW=1:Level, standX=TRUE)
  preds <- predict(KDSNpart, newx=as.matrix(trees[-treeTestInd[[i]], -3]))
  errKDSN[i] <- sqrt(mean((trees[-treeTestInd[[i]], "Volume"]-preds)^2))
}
errKDSN
mean((errKDSN - err))

# Fine tuning
Level <- 3
for(i in 1:10) {
  KDSNpart <- fitKDSN(y=trees[treeTestInd[[i]], "Volume"], X=as.matrix(trees[treeTestInd[[i]], -3]), 
                      levels=Level, Dim=round(seq(dim(trees)[1], sqrt(dim(trees)[1]), length.out=Level)), 
                      sigma=seq(quantEuklid[1], quantEuklid[2], length.out=Level), lambda=seq(10^-1, 10^-10, length.out=Level), 
                      alpha=rep(0, Level), 
                      info=FALSE, seedW=1:Level, standX=TRUE)
  preds <- predict(KDSNpart, newx=as.matrix(trees[-treeTestInd[[i]], -3]))
  errKDSN[i] <- sqrt(mean((trees[-treeTestInd[[i]], "Volume"]-preds)^2))
}
mean(errKDSN)

errKDSN1 <- vector("numeric", 10)
for(i in 1:10) {
  KDSNpart <- fineTuneKDSN(KDSNpart, 
                           y=matrix(trees[treeTestInd[[i]], "Volume"], ncol=1), 
                           X=as.matrix(trees[treeTestInd[[i]], -3]),
                           fineTuneIt=100, info=FALSE, seedInit = 0)
  preds <- predict(KDSNpart, newx=as.matrix(trees[-treeTestInd[[i]], -3]))
  errKDSN1[i] <- sqrt(mean((trees[-treeTestInd[[i]], "Volume"]-preds)^2))
  cat("iter = ", i, "\n")
}
errKDSN1
stopifnot(mean(errKDSN1) < mean(errKDSN))
