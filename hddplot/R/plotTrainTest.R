"plotTrainTest" <-
function(x, nfeatures, cl, traintest,
           titles=c("A: I/II (train with I, scores are for II)",
             "B: II/I (train with II, scores are for I)")){
    oldpar <- par(mfrow=c(1,2), pty="s")
    on.exit(par(oldpar))
    if(length(nfeatures)==1)nfeatures <- rep(nfeatures,2)
    traintest <- factor(traintest)
    train <- traintest==levels(traintest)[1]
    testset <- traintest==levels(traintest)[2]
    cl1 <- cl[train]
    cl2 <- cl[testset]
    nf1 <- nfeatures[1]
    ord1 <- orderFeatures(x, cl, subset=train)
    df1 <- data.frame(t(x[ord1[1:nf1], train]))
    df2 <- data.frame(t(x[ord1[1:nf1], testset]))
    df1.lda <- lda(df1, cl1)
    scores <- predict(df1.lda, newdata=df2)$x
    scoreplot(scorelist=list(scores=scores, cl=cl2,
             nfeatures=nfeatures[1], other=NULL, cl.other=NULL),
           prefix.title="")
    mtext(side=3, line=2, titles[1], adj=0)
    nf2 <- nfeatures[2]
    ord2 <- orderFeatures(x, cl, subset=testset)
    df2 <- data.frame(t(x[ord2[1:nf2], testset]))
    df1 <- data.frame(t(x[ord2[1:nf2], train]))
    df2.lda <- lda(df2, cl2)
    scores <- predict(df2.lda, newdata=df1)$x
    scoreplot(scorelist=list(scores=scores, cl=cl1,
             nfeatures=nfeatures[2], other=NULL, cl.other=NULL),
           prefix.title="")
    mtext(side=3, line=2, titles[2], adj=0)
  }

