## This chapter uses data and functions from some packages that are
## not automatically installed when installing
## ChemometricsWithR. The script checks their presence and in case they
## are absent does not execute the corresponding code.
if (!require("mclust")) {
  mclust.present <- FALSE
  cat("Package mclust not available - some code may not run.\nInstall it by typing 'install.packages(\"mclust\")'")
} else {
  mclust.present <- TRUE
}
if (!require("sfsmisc")) {
  sfsmisc.present <- FALSE
  cat("Package sfsmisc not available - some code may not run.\nInstall it by typing 'install.packages(\"sfsmisc\")'")
} else {
  sfsmisc.present <- TRUE
}
if (!require("class")) {
  class.present <- FALSE
  cat("Package class not available - some code may not run.\nInstall it by typing 'install.packages(\"class\")'")
} else {
  class.present <- TRUE
}
if (!require("rda")) {
  rda.present <- FALSE
  cat("Package rda not available - some code may not run.\nInstall it by typing 'install.packages(\"rda\")'")
} else {
  rda.present <- TRUE
}
if (!require("rpart")) {
  rpart.present <- FALSE
  cat("Package rpart not available - some code may not run.\nInstall it by typing 'install.packages(\"rpart\")'")
} else {
  rpart.present <- TRUE
}
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

data(wines, package = "ChemometricsWithRData")
odd <- seq(1, nrow(wines), by = 2)
even <- seq(2, nrow(wines), by = 2)
wines.trn <- wines[odd,]
wines.tst <- wines[even,]
## discriminant analysis
wines.counts <- table(vintages[odd])
ngroups <- length(wines.counts)
wines.groups <- split(as.data.frame(wines.trn), 
                      vintages[odd])
wines.covmats <- lapply(wines.groups, cov)
wines.wcovmats <- lapply(1:ngroups,
                         function(i, x, y) x[[i]]*y[i],
                         wines.covmats, wines.counts)
wines.pooledcov <- Reduce("+", wines.wcovmats) /
  (nrow(wines.trn) - ngroups)
wines.pooledcov2 <- matrix(0, ncol(wines), ncol(wines))
for (i in 1:3) {
  wines.pooledcov2 <- wines.pooledcov2 +
    cov(wines.groups[[i]]) * nrow(wines.groups[[i]])
}
wines.pooledcov2 <- 
  wines.pooledcov2 / (nrow(wines.trn) - ngroups)
range(wines.pooledcov2 - wines.pooledcov)

distances <- 
  sapply(1:ngroups, 
         function(i, samples, means, covs)
         mahalanobis(samples, colMeans(means[[i]]), covs),
         wines.trn, wines.groups, wines.pooledcov)
trn.pred <- apply(distances, 1, which.min)
table(vintages[odd], trn.pred)

distances <- 
  sapply(1:ngroups,
         function(i, samples, means, covs)
         mahalanobis(samples, colMeans(means[[i]]), covs),
         wines.tst, wines.groups, wines.pooledcov)
tst.pred <- apply(distances, 1, which.min)
table(vintages[even], tst.pred)

wines.ldamod <- lda(wines.trn, grouping = vintages[odd],
                    prior = rep(1,3)/3)
wines.lda.testpred <- predict(wines.ldamod, new = wines.tst)
table(vintages[even], wines.lda.testpred$class)

plot(wines.ldamod, xlim = c(-7, 6), col = as.integer(vintages[odd]))

wines.ldamod <- lda(wines.trn, grouping = vintages[odd],
                    prior = rep(1,3)/3, CV = TRUE)
table(vintages[odd],wines.ldamod$class)

X <- wines[vintages != "Barolo", c(7, 13)]
vint <- factor(vintages[vintages != "Barolo"])
wines.counts <- table(vint)
wines.groups <- split(as.data.frame(X), vint)
WSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x) {
                       crossprod(scale(x, scale = FALSE))}))
BSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x, y) {
                       nrow(x) * tcrossprod(colMeans(x) - y)},
                     colMeans(X)))
FLDA <- eigen(solve(WSS, BSS))$vectors[,1]
FLDA / FLDA[1]

wines.covmats <- lapply(wines.groups, cov)
wines.wcovmats <- lapply(1:length(wines.groups),
                         function(i, x, y) x[[i]]*y[i],
                         wines.covmats, wines.counts)
wines.pcov12 <- Reduce("+", wines.wcovmats) / (length(vint) - 2)
MLLDA <-
  solve(wines.pcov12,
        apply(sapply(wines.groups, colMeans), 1, diff))
MLLDA / MLLDA[1]

## Three-group case, not explicitly written in the book
wines.trn <- wines[odd, c(7, 13)]
wines.tst <- wines[even, c(7, 13)]
wines.counts <- table(vintages[odd])
wines.groups <- split(as.data.frame(wines.trn), vintages[odd])
WSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x) {
                       crossprod(scale(x, scale = FALSE))}))
BSS <- 
  Reduce("+", lapply(wines.groups,
                     function(x, y) {
                       nrow(x) * tcrossprod(colMeans(x) - y)},
                     colMeans(X)))
FLDA <- eigen(solve(WSS, BSS))$vectors[,1]

x <- seq(.4, 5.4, length = 251)
y <- seq(250, 1750, length = 251)
gridXY <- cbind(rep(x, each = length(y)), rep(y, length(x)))

scores <- gridXY %*% FLDA
meanscores <- t(sapply(wines.groups, colMeans)) %*% FLDA
Fdistance <- outer(c(scores), c(meanscores), 
                   FUN = function(x, y) abs(x - y))
Fclassif <- apply(Fdistance, 1, which.min)

contour(x, y,
        matrix(Fclassif, nrow = length(x), byrow = TRUE),
        main = "Fisher LDA", drawlabels = FALSE,
        xlab = "flavonoids", ylab = "proline")
points(wines.tst, col = as.integer(vintages[even]),
       pch = as.integer(vintages[even]))

wines.ldamod <- lda(wines.trn, 
                    grouping = vintages[odd],
                    prior = rep(1,3)/3)
lda.2Dclassif <- predict(wines.ldamod, newdata = gridXY)$class
contour(x, y,
        matrix(as.integer(lda.2Dclassif), 
               nrow = length(x), byrow = TRUE),
        main = "LDA", drawlabels = FALSE,
        xlab = "flavonoids", ylab = "proline")
points(wines.tst, col = as.integer(vintages[even]),
       pch = as.integer(vintages[even]))

## QDA
wines.trn <- wines[odd, c(7, 13)]
wines.tst <- wines[even,  c(7, 13)]
wines.groups <- split(as.data.frame(wines.trn), vintages[odd])
wines.covmats <- lapply(wines.groups, cov)
ngroups <- length(wines.groups)
distances <- sapply(1:ngroups,
                    function(i, samples, means, covs) {
                      mahalanobis(samples,
                                  colMeans(means[[i]]),
                                  covs[[i]]) },
                    wines.tst, wines.groups, wines.covmats)
test.pred <- apply(distances, 1, which.min)
table(vintages[even], test.pred)
qda.mahal.dists <- sapply(1:ngroups,
                          function(i, samples, means, covs) {
                            mahalanobis(samples,
                                        colMeans(means[[i]]),
                                        covs[[i]]) },
                          gridXY, wines.groups, wines.covmats)
qda.2Dclassif <- apply(qda.mahal.dists, 1, which.min)
contour(x, y, 
        matrix(qda.2Dclassif, nrow = length(x), byrow = TRUE),
        main = "QDA", drawlabels = FALSE,
        xlab = "flavonoids", ylab = "proline")
points(wines.tst, col = as.integer(vintages[even]),
       pch = as.integer(vintages[even]))

wines.qda <- qda(wines[odd,], vintages[odd], 
                 prior = rep(1, 3)/3)
test.qdapred <- predict(wines.qda, newdata = wines[even,])
table(vintages[even], test.qdapred$class)


if (mclust.present) {
  ## MBDA
  wines.mclustDA <- mclustDA(train = list(data = wines[odd,],
                               labels = vintages[odd]),
                             test = list(data = wines[even,],
                               labels = vintages[even]),
                             G = 1:5)
  wines.mclustDA
  
  wines.mclust2D <- mclustDAtrain(wines[odd, c(7, 13)],
                                  vintages[odd], G = 1:5)
  wines.mclust2Dpred <- mclustDAtest(gridXY, wines.mclust2D)  
  contour(x, y,
          matrix(apply(wines.mclust2Dpred, 1, which.max),
                 nrow = length(x), byrow = TRUE),
          main = "MBDA", drawlabels = FALSE,
          xlab = "flavonoids", ylab = "proline")
  points(wines[even, c(7, 13)],
         col = as.integer(vintages[even]),
         pch = as.integer(vintages[even]))
}

if (sfsmisc.present & rda.present) {
  ## DDA
  data(prostate, package = "ChemometricsWithRData")
  prost <- prostate[prostate.type != "bph", 1:1000]
  prost.type <- factor(prostate.type[prostate.type != "bph"])
  odd <- seq(1, length(prost.type), by = 2)
  even <- seq(2, length(prost.type), by = 2)
  prost.dlda <-
    dDA(prost[odd,], as.integer(prost.type)[odd], pool = TRUE)
  prost.dldapred <- predict(prost.dlda, prost[even,])
  table(prost.type[even], prost.dldapred)
  
  prost.dqda <-
    dDA(prost[odd,], as.integer(prost.type)[odd], pool = FALSE)
  prost.dqdapred <- predict(prost.dqda, prost[even,])
  table(prost.type[even], prost.dqdapred)

  ## RDA
  prost.rda <-
    rda(t(prost[odd,]), as.integer(prost.type)[odd],
        delta = seq(0, .4, length = 5),
        alpha = seq(0, .4, length = 5))
  prost.rda
  
  prost.rdacv <-                                                               
    rda.cv(prost.rda, t(prost[odd,]),                                          
           as.integer(prost.type)[odd])
  prost.rdapred <-
    predict(prost.rda,
            t(prost[odd,]), as.integer(prost.type)[odd],
            t(prost[even,]), alpha = .2, delta = 0)
  table(prost.type[even],prost.rdapred)
}

if (class.present & e1071.present) {
  ## Nearest neighbours
  odd <- seq(1, nrow(wines), by = 2)
  even <- seq(2, nrow(wines), by = 2)
  wines.sc <- scale(wines, scale = apply(wines[odd,], 2, sd), 
                    center = colMeans(wines[odd,]))
  dist2sample2a <- mahalanobis(wines.sc[odd,], wines.sc[2,], 
                               diag(13))
  dist2sample2b <- mahalanobis(wines[odd,], wines[2,],
                               diag(apply(wines[odd,], 2, sd)^2))
  range(dist2sample2a - dist2sample2b)
  
  nearest.classes <- vintages[odd][order(dist2sample2a)]
  nearest.classes[1:10]
  
  dist2sample2 <- mahalanobis(wines[odd,], wines[2,],
                              cov(wines[odd,]))
  nearest.classes <- vintages[odd][order(dist2sample2)]
  nearest.classes[1:10]
  
  X <- scale(wines, scale = apply(wines[odd,], 2, sd),
             center = colMeans(wines[odd,]))
  knn(X[odd,], X[68,], cl = vintages[odd], k = 4)
  knn(X[odd,], X[68,], cl = vintages[odd], k = 4)
  
  knn(X[odd,], X[68,], cl = vintages[odd], k = 4, l = 3)
  
  wines.knnresult <- rep(0, 10)
  for (i in 1:10) {
    wines.knncv <- knn.cv(X[odd,], vintages[odd], k = i)
    wines.knnresult[i] <-
      sum(diag(table(vintages[odd], wines.knncv))) }
  100 * wines.knnresult / length(odd)
  
  set.seed(7)
  knn.tuned <- tune.knn(X[odd,], vintages[odd], k = 1:10)
  knn.tuned
  
  plot(knn.tuned)
  ## Next bit takes some time...
  ## bestKs <- rep(0, 1000)
  ## for (i in 1:1000)
  ##   bestKs[i] <- tune.knn(X1, vintages[odd], 
  ##                         k = 1:10)$best.parameters[1,1]
  ## hist(bestKs)
}

if (rpart.present) {
  ## trees - rpart
  wines.df <- data.frame(vint = vintages, wines[,c(7, 13)])
  wines.rpart <- rpart(vint ~ ., subset = odd,
                       data = wines.df, method = "class")
  wines.rpart
  
  plot(wines.rpart, margin = .12)
  text(wines.rpart, use.n = TRUE)
  plot(wines[even, c(7, 13)], pch = as.integer(vintages[even]),
       col = as.integer(vintages[even]))
  segments(wines.rpart$splits[1,4], par("usr")[3],
           wines.rpart$splits[1,4], par("usr")[4], lty = 2)
  segments(wines.rpart$splits[1,4], wines.rpart$splits[2,4],
           par("usr")[2], wines.rpart$splits[2,4], lty = 2)
  
  wines.df <- data.frame(wines, vint = vintages)
  wines.rpart <- rpart(vint ~ ., subset = odd,
                       data = wines.df, method = "class")
  plot(wines.rpart, margin = .1)
  text(wines.rpart, use.n = TRUE)
  
  wines.rpart.predict <- predict(wines.rpart, 
                                 newdata = wines.df[even,])
  wines.rpart.predict[31:34,]
  
  matplot(wines.rpart.predict)
  
  table(vintages[even], 
        predict(wines.rpart, newdata = wines.df[even,], 
                type = "class"))
  
  X <- wines[odd,c(7, 13)]
  Ginis <- matrix(0, nrow(X), 2)
  splits.flav <- sort(X[,1])
  splits.prol <- sort(X[,2])
  for (i in 1:nrow(X)) {
    Ginis[i,1] <- gini(X[,1], vintages[odd], splits.flav[i])
    Ginis[i,2] <- gini(X[,2], vintages[odd], splits.prol[i])
  }
  matplot(Ginis, type = "l", lty = 1:2, col = 2:1)
  legend("top", legend = c("flavonoids", "proline"),
         col = 2:1, lty = 1:2)
  
  apply(Ginis, 2, which.min)
  
  odd <-   odd <- seq(1, length(prost.type), by = 2)
  even <- seq(2, length(prost.type), by = 2)
  prost.df <- data.frame(type = prost.type, prost = prost)
  prost.rprt <- 
    rpart(type ~ ., data = prost.df, subset = odd,
          control = rpart.control(cp = 0, minsplit = 0))
  prost.rprtpred <- 
    predict(prost.rprt, newdata = prost.df[even,])
  table(prost.type[even], classmat2classvec(prost.rprtpred))
  
  printcp(prost.rprt)
  
  plotcp(prost.rprt)
  
  prost.rprt2 <- 
    rpart(type ~ ., data = prost.df, subset = odd,
          control = rpart.control(cp = 0.12))
  prost.rprt2pred <- 
    predict(prost.rprt2, newdata = prost.df[even,])
  table(prost.type[even], classmat2classvec(prost.rprt2pred))
}

if (e1071.present) {
  ## SVMs
  wines.df <- data.frame(vint = factor(vintages[vintages != "Barolo"]), 
                         wines[vintages != "Barolo",])
  odd <- seq(1, nrow(wines.df), by = 2)
  even <- seq(2, nrow(wines.df), by = 2)
  wines.svm <- svm(vint ~ ., data = wines.df, subset = odd)
  wines.svmpred <- predict(wines.svm, newdata = wines.df[even,])
  table(vint[even], wines.svmpred)
  
  odd <-   odd <- seq(1, length(prost.type), by = 2)
  even <- seq(2, length(prost.type), by = 2)
  prost <- prostate[prostate.type != "bph",1:1000]
  prost.type <- factor(prostate.type[prostate.type != "bph"])
  prost.df <- data.frame(type = prost.type, prost = prost)
  prost.svm <- svm(type ~ ., data = prost.df, subset = odd,
                   cross = 10)
  summary(prost.svm)
  
  prost.svmpred <- predict(prost.svm, newdata = prost.df[even,])
  table(prost.type[even], prost.svmpred)
  
  ## Multiclass SVMs
  odd <- seq(1, nrow(wines), by = 2)
  even <- seq(2, nrow(wines), by = 2)
  wines.df <- data.frame(vintages = vintages, wines[,c(7,13)])
  wines.svm <- svm(vintages ~ ., data = wines.df, subset = odd)
  
  ## The optimisation shows variable results depending on the state of
  ## the random generator...
  set.seed(13)
  wines.dfodd <- wines.df[odd,]
  wines.bestsvm <-
    best.svm(vintages ~ ., data = wines.dfodd,
             kernel = "polynomial", 
             coef0 = seq(-.5, .5, by = .1),
             gamma = 2^(-1:1), cost = 2^(2:4))
  wines.bestsvmpred <- 
    predict(wines.bestsvm, newdata = wines.df[even,])
  sum(wines.bestsvmpred == vintages[even])
  
  plot(wines.svm, wines.df[odd,], proline ~ flavonoids, grid = 500)
  plot(wines.bestsvm, wines.df[odd,], proline ~ flavonoids, grid = 500)
}

if (nnet.present & e1071.present) {
  ## ANNs
  X <- scale(wines, scale = apply(wines[odd,], 2, sd),
             center = colMeans(wines[odd,]))
  w.df <- data.frame(vintage = vintages, wines = X)
  set.seed(7)
  wines.nnet <- nnet(vintage ~ ., data = w.df,
                     size = 4, subset = odd)
  
  training.pred <- predict(wines.nnet, type = "class")
  sum(diag(table(vintages[odd],
                 training.pred))) / length(odd)
  
  table(vintages[even],
        predict(wines.nnet, w.df[even,], type = "class"))
  
  wines.nnetmodels <- 
    tune.nnet(vintage ~ ., data = w.df[odd,], size = 1:8)
  summary(wines.nnetmodels)
  
  set.seed(7)
  best.wines.nnet <- best.nnet(vintage ~ ., data = w.df[odd,],
                               size = 1:8)
  table(vintages[even],
        predict(best.wines.nnet, w.df[even,], type = "class"))
}
