"compareTreecalcs" <-
function (x = yesno ~ ., data = DAAG::spam7, cp = 0.00025, fun = c("rpart",
                                                         "randomForest"))
{
  rpart.out <- try(requireNamespace("rpart"), silent = TRUE)
  rf.out <- try(requireNamespace("randomForest"), silent = TRUE)
  rpartCheck <- !is.logical(rpart.out) | (rpart.out == FALSE)
  rfCheck <- !is.logical(rf.out) | (rf.out == FALSE)
  if(rpartCheck)print("Error: package rpart is not installed properly")
  if(rfCheck)print("Error: package randomForest is not installed properly")
  if(rpartCheck | rfCheck) return()
  yvar <- all.vars(x)[1]
  m <- dim(data)[1]
  train <- sample((1:m), m%/%2)
  dftrain <- data[train, ]
  dftest <- data[-train, ]
  if ("rpart" %in% fun) {
    df.rpart <- rpart::rpart(x, data = dftrain, cp = cp)
    cptable <- df.rpart$cptable
    err.root <- df.rpart$frame$dev[1]/df.rpart$frame$n[1]
    if("xerror" %in% colnames(cptable)){
       xerror <- cptable[, "xerror"]
       xstd <- cptable[, "xstd"]
       CP <- cptable[, "CP"]
       if (which.min(xerror)==length(xerror))
         write.table(data.frame(c("Minimum cv error was for largest tree",
                       "Perhaps try a smaller value of cp")),
                     row.names=FALSE, col.names=FALSE, quote=FALSE)
     }
       else {
       write.table(data.frame("Warning: Calculations ceased at the root node"),
                      row.names=FALSE, col.names=FALSE, quote=FALSE)
       xerror <- numeric(0)
     }
    nREmin <- max(as.numeric(which.min(xerror)),0)
    if (length(nREmin)>0 & nREmin != length(xerror)) {
      re.min <- xerror[nREmin]
      cv.min <- as.vector(re.min * err.root)
      selim <- min(xerror + xstd)
      nSErule <- min(seq(along = xerror)[xerror <= selim])
      cv.selim <- as.vector(err.root * xerror[nSErule])
      cp.selim <- mean(CP[(nSErule - 1):nSErule])
      cp.remin <- mean(CP[(nREmin - 1):nREmin])
      df.rpart <- rpart::prune(df.rpart, cp = cp.remin)
      hat <- predict(df.rpart, newdata = dftest, type = "class")
      tab <- table(hat, dftest[, yvar])
      testerr.estmin <- sum(tab[row(tab) != col(tab)])/sum(tab)
      df.rpart <- rpart::prune(df.rpart, cp = cp.selim)
      hat <- predict(df.rpart, newdata = dftest, type = "class")
      tab <- table(hat, dftest[, yvar])
      testerr.SE <- sum(tab[row(tab) != col(tab)])/sum(tab)
    }
    else {
      cv.selim <- NA
      cv.min <- NA
      testerr.SE <- NA
      testerr.estmin <- NA
      nSErule <- NA
      nREmin <- NA
    }
  }
  else {
    cv.selim <- NULL
    cv.min <- NULL
    testerr.SE <- NULL
    testerr.estmin <- NULL
    nSErule <- NULL
    nREmin <- NULL
  }
  if ("randomForest" %in% fun) {
    y <- dftrain[, yvar]
    ynum <- match(yvar, names(dftrain))
    df.rf <- randomForest::randomForest(x = dftrain[, -ynum], y = y)
    hat.dfcv <- predict(df.rf, type = "response")
    tab <- table(hat.dfcv, dftrain[, yvar])
    rfcvA <- sum(tab[row(tab) != col(tab)])/sum(tab)
    hat.dftest <- predict(df.rf, newdata = dftest[, -ynum],
                          type = "response")
    tab <- table(hat.dftest, dftest[, yvar])
    rftest <- sum(tab[row(tab) != col(tab)])/sum(tab)
  }
  else {
    rfcvA <- NULL
    rftest <- NULL
  }
  c(rpSEcvI = cv.selim, rpcvI = cv.min, rpSEtest = testerr.SE,
    rptest = testerr.estmin, nSErule = nSErule, nREmin = nREmin,
    rfcvI = rfcvA, rftest = rftest)
}

