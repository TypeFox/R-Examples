RPChooseSS <-
    function(XTrain #n by p training data matrix
           , YTrain #n vector of the classes of the training samples
           , XVal #n.val by p validation data matrix
           , YVal #n.val vector of the classes of the validation samples
           , XTest  #n.test by p test data matrix
           , d      #dimension to project into
           , B2 = 100 #block size
           , base = "LDA" # base classifier, eg "knn","LDA","QDA" or other
           , k = c(3,5) # possible k if base = "knn"
           , projmethod = "Haar" # projection distribution eg "Haar", "axis"
         ,   ... )
        {
    n <- length(YTrain)
    p <- ncol(XTrain)
    n.val <- length(YVal)
    w <- n.val
    if (base == "LDA") {
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            weight.test <- mean(predict(lda(x = XTrain%*%RP, grouping = YTrain), XVal%*%RP)$class != YVal, na.rm = TRUE)
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
            }
        }
        Val.Class <- as.numeric(predict(lda(x = XTrain%*%RP1, grouping = YTrain), XVal%*%RP1)$class)
        Test.Class <- as.numeric(predict(lda(x = XTrain%*%RP1, grouping = YTrain), XTest%*%RP1)$class)
    }
    if (base == "QDA") {
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            weight.test <- mean(predict(qda(x = XTrain%*%RP, grouping = YTrain), XVal%*%RP)$class != YVal, na.rm = TRUE)
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
            }
        }
        Val.Class <- as.numeric(predict(qda(x = XTrain%*%RP1, grouping = YTrain), XVal%*%RP1)$class)
        Test.Class <- as.numeric(predict(qda(x = XTrain%*%RP1,  grouping = YTrain), XTest%*%RP1)$class)
    }
    if (base == "knn"){
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            kcv.voteRP <- sapply(k, function(x) {mean(knn(XTrain%*%RP, XVal%*%RP, YTrain, x) != YVal, na.rm = TRUE)})
            weight.test <- min(kcv.voteRP)
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
                k1 <- order(kcv.voteRP)[1]
            }
        }
        Val.Class <- as.numeric(knn(XTrain%*%RP1, XVal%*%RP1, YTrain, k1))
        Test.Class <- as.numeric(knn(XTrain%*%RP1, XTest%*%RP1, YTrain, k1))
    }
    if (base == "Other") {
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            weight.test <- mean(Other.classifier(x = XTrain%*%RP, grouping = YTrain, XVal%*%RP)$class != YVal, na.rm = TRUE)
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
            }
        }
        Val.Class <- as.numeric(Other.classifier(x = XTrain%*%RP1, grouping = YTrain, XVal%*%RP1)$class)
        Test.Class <- as.numeric(Other.classifier(x = XTrain%*%RP1, grouping = YTrain, XTest%*%RP1)$class)
    }
    return(c(Val.Class, Test.Class)) 
}
