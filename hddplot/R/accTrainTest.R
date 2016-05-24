"accTrainTest" <-
function(x=matrix(rnorm(1000), ncol=20), cl=factor(rep(1:3,c(7,9,4))),
         traintest=divideUp(cl, nset=2), nfeatures=NULL, print.acc=FALSE,
         print.progress=TRUE){
    traintest <- factor(traintest)
    train <- traintest==levels(traintest)[1]
    testset <- traintest==levels(traintest)[2]
    cl1 <- cl[train]
    cl2 <- cl[testset]
    ng1 <- length(cl1)
    ng2 <- length(cl2)
    maxg <- max(c(ng1-length(unique(cl1))-2,
                  ng2-length(unique(cl2))-2))
    if(is.null(nfeatures)){
      max.features <- maxg
      nfeatures <- 1:max.features
    } else
    {
      if(max(nfeatures)>maxg)nfeatures <- nfeatures[nfeatures<=maxg]
      max.features <- max(nfeatures)
    }
    ord1 <- orderFeatures(x, cl, subset=train)[1:max.features]
    ord2 <- orderFeatures(x, cl, subset=testset)[1:max.features]
    ord <- unique(c(ord1, ord2))
    sub1 <- match(ord1, ord)
    sub2 <- match(ord2, ord)
    df1 <- data.frame(t(x[ord, train]))
    df2 <- data.frame(t(x[ord, testset]))
    acc1 <- acc2 <- numeric(max(nfeatures))
    for(i in nfeatures){
      if(print.progress)cat(paste(i, ":", sep=""))
      df1.lda <- lda(df1[, sub1[1:i], drop=FALSE], cl1)
      hat2 <- predict(df1.lda, newdata=df2[, sub1[1:i], drop=FALSE])$class
      tab <- table(hat2, cl2)
      acc1[i] <- sum(tab[row(tab)==col(tab)])/sum(tab)
      df2.lda <- lda(df2[, sub2[1:i], drop=FALSE], cl2)
      hat1 <- predict(df2.lda, newdata=df1[, sub2[1:i], drop=FALSE])$class
      tab <- table(hat1, cl1)
      acc2[i] <- sum(tab[row(tab)==col(tab)])/sum(tab)
    }
    cat("\n")
    if(print.acc){
      print(round(acc1,2))
      print(round(acc2,2))
    }
    maxacc1 <- max(acc1)
    maxacc2 <- max(acc2)
    sub1 <- match(maxacc1, acc1)
    sub2 <- match(maxacc2, acc2)
    nextacc1 <- max(acc1[acc1<1])
    nextacc2 <- max(acc1[acc1<2])
    lower1 <- maxacc1-sqrt(nextacc1*(1-nextacc1)/ng1)
    lower2 <- maxacc2-sqrt(nextacc2*(1-nextacc2)/ng2)
    lsub1 <- min((1:ng1)[acc1>lower1])
    lsub2 <- min((1:ng2)[acc2>lower2])
    lower <- c("Best accuracy, less 1SD  ",
               paste(paste(round(c(lower1, lower2),2), c(lsub1, lsub2),
                           sep=" ("), " features)   ", sep=""))
    best <- c("Best accuracy",
              paste(paste(round(c(maxacc1, maxacc2),2), c(sub1, sub2),
                          sep=" ("), " features)", sep=""))
    acc.df <- cbind(lower, best)
    dimnames(acc.df) <- list(c("Training/test split",
                               "I (training) / II (test)    ",
                               "II (training) / I (test)    "),c("",""))
    print(acc.df, quote=FALSE)
    invisible(list(sub1.2=ord1, acc1.2=acc1, sub2.1=ord2, acc2.1=acc2))
  }

