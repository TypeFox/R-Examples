### R code from vignette source 'bmrm.Rnw'

###################################################
### code chunk number 1: bmrm.Rnw:27-34
###################################################
  require(bmrm)
  # -- Create a 2D dataset with the first 2 features of iris, and binary labels
  x <- data.matrix(iris[1:2])
  y <- c(-1,1,1)[iris$Species]

  # -- Add a constant dimension to the dataset to learn the intercept
  x <- cbind(x,1)


###################################################
### code chunk number 2: bmrm.Rnw:37-51
###################################################
  train.prediction.model <- function(...) {
    m <- bmrm(...)
    m$f <- x %*% m$w
    m$y <- sign(m$f)
    m$contingencyTable <- table(y,m$y)
    return(m)
  }

  # -- train models with maxMarginLoss and fbetaLoss
  models <- list(
    svm_L1 = train.prediction.model(hingeLoss(x,y),LAMBDA=0.01,regfun='l1'),
    svm_L2 = train.prediction.model(hingeLoss(x,y),LAMBDA=0.1,regfun='l2'),
    f1_L1 = train.prediction.model(fbetaLoss(x,y),LAMBDA=0.01,regfun='l1')
  )


###################################################
### code chunk number 3: bmrm.Rnw:56-72
###################################################
  # -- Plot the dataset and the predictions
  layout(matrix(1:2,1,2))
  plot(x,pch=20+y,main="dataset & hyperplanes")
  legend('bottomright',legend=names(models),col=seq_along(models),lty=1,cex=0.75,lwd=3)
  for(i in seq_along(models)) {
    m <- models[[i]]
    if (m$w[2]!=0) abline(-m$w[3]/m$w[2],-m$w[1]/m$w[2],col=i,lwd=3)
  }

  rx <- range(na.rm=TRUE,1,unlist(lapply(models,function(e) nrow(e$log))))
  ry <- range(na.rm=TRUE,0,unlist(lapply(models,function(e) e$log$epsilon)))
  plot(rx,ry,type="n",ylab="epsilon gap",xlab="iteration",main="evolution of the epsilon gap")
  for(i in seq_along(models)) {
    m <- models[[i]]
    lines(m$log$epsilon,type="o",col=i,lwd=3)
  }


