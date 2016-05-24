.packageName <- 'upclass'

noupclassifymodel <- function (Xtrain, cltrain, Xtest, cltest = NULL, modelName = "EEE", ...) 
{
  functioncall <- match.call(expand.dots = TRUE)
  if (is.matrix(Xtrain)) {
    Ntrain <- nrow(Xtrain)
  }
  else {
    Ntrain <- length(Xtrain)
  }
  if (is.matrix(Xtest)) {
    Ntest <- nrow(Xtest)
  }
  else {
    Ntest <- length(Xtest)
  }
  ltrain <- unmap(cltrain)
  grps<-names(table(cltrain)) # need this for the unmap command
  G <- ncol(ltrain)
  if (!is.null(cltest)) {
    ltest <-unmap(cltest,groups=grps)   
  }
  if (is.matrix(Xtrain) & (!is.matrix(Xtest))) {
    Xtest <- matrix(Xtest, 1, length(Xtest))
    if (!is.null(cltest)) {
      ltest <- matrix(ltest, 1, length(ltest))
    }
  }
  if (is.matrix(Xtrain)) {
    d <- ncol(Xtrain)
  }
  else {
    d <- 1
  }
  G <- dim(ltrain)[2]
  fitm <- mstep(modelName = modelName, data = Xtrain, z = ltrain)
  fite <- do.call("estep", c(list(data = Xtest), fitm))
  z <- fite$z
  mTau.train <- matrix(log(fitm$parameters$pro), Ntrain, fitm$G, 
                       byrow = TRUE)
  lDensity.train <- do.call("cdens", c(list(data = Xtrain, 
                                            logarithm = TRUE), fitm))
  sum.train <- mTau.train + lDensity.train
  mat.train <- ltrain * sum.train
  ll <- sum(mat.train)
  res <- list()
  res$call <- functioncall
  res$Ntrain <- Ntrain
  res$Ntest <- Ntest
  res$d <- d
  res$G <- G
  res$modelName <- modelName
  res$parameters <- fitm$parameters
  
    fitetrain <- do.call("estep", c(list(data = Xtrain), 
                                    fitm))
    ztrain <- fitetrain$z
    cltrain <- map(ztrain)
    tab0 <- table(map(ltrain), cltrain)
    tab <- matrix(0, G, G)
    rownames(tab) <- colnames(tab) <- 1:G
    tab[rownames(tab0), colnames(tab0)] <- tab0
    rate <- sum(tab) - sum(diag(tab))
    Brier <- sum((ltrain - ztrain)^2)*100/(2*Ntrain)
    res$train <- list()
    res$train$z <- ztrain
    res$train$cl <- cltrain
    res$train$misclass <- rate
    res$train$rate <-rate*100/length(cltrain)
    res$train$Brierscore <- Brier
    res$train$tab <- tab
    cl <- map(z)
    res$test <- list()
    res$test$z <- z
    res$test$cl <- cl
    if (!is.null(cltest)) {
      cltest <- map(ltest)
      tab0 <- table(cltest, cl)
      tab <- matrix(0, G, G)
      rownames(tab) <- colnames(tab) <- 1:G
      tab[rownames(tab0), colnames(tab0)] <- tab0
      rate <- sum(tab) - sum(diag(tab))
      Brier <- sum((ltest - z)^2)*100/(2*Ntest)
      res$test$misclass <- rate
      res$test$rate <-rate*100/length(cltest)
      res$test$Brierscore <- Brier
      res$test$tab <- tab
    }
  bic.all <- bic(modelName, ll, fitm$n, fitm$d, fitm$G)
  res$ll <- ll
  res$bic <- bic.all
  res
}