RPChoose <-
    function(XTrain #n by p trining data matrix
           , YTrain #n vector of the classes of the trining samples
           , XTest  #n.test by p test data matrix
           , d      #dimension to project into
           , B2 = 100 #block size
           , base = "LDA" # base classifier, eg "knn","LDA","QDA" or other
           , k = c(3,5) # possible k if base = "knn"
           , projmethod = "Haar" # projection distribution eg "Haar", "axis"
           , estmethod = "resub"
         ,   ... )
        {
    n <- length(YTrain)
    p <- ncol(XTrain)
    w <- n
    if (base == "LDA") {
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            if(estmethod == "resub"){
                weight.test <- mean(predict(lda(x = XTrain%*%RP, grouping = YTrain), XTrain%*%RP)$class != YTrain, na.rm = TRUE)
            }
            if(estmethod == "loo"){
                weight.test <- mean(lda(x = XTrain%*%RP, grouping = YTrain, CV = TRUE)$class != YTrain, na.rm = TRUE)
            }
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
            }
        }
        if(estmethod == "resub"){Train.Class <- as.numeric(predict(lda(x = XTrain%*%RP1, grouping = YTrain), XTrain%*%RP1)$class)}
        if(estmethod == "loo"){Train.Class <- as.numeric(lda(x = XTrain%*%RP1, grouping = YTrain, CV = TRUE)$class)}
        Test.Class <- as.numeric(predict(lda(x = XTrain%*%RP1, grouping = YTrain), XTest%*%RP1)$class)
    }
    if (base == "QDA") {
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            if(estmethod == "resub"){
                weight.test <- mean(predict(qda(x = XTrain%*%RP, grouping = YTrain), XTrain%*%RP)$class != YTrain, na.rm = TRUE)
            }
            if(estmethod == "loo"){
                weight.test <- mean(qda(x = XTrain%*%RP, grouping = YTrain,CV = TRUE)$class != YTrain, na.rm = TRUE)
            }
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
            }
        }
        if(estmethod == "resub"){Train.Class <- as.numeric(predict(qda(x = XTrain%*%RP1, grouping = YTrain), XTrain%*% RP1)$class)}
        if(estmethod == "loo"){Train.Class <- as.numeric(qda(x = XTrain%*%RP1, grouping = YTrain, CV = TRUE)$class)}
        Test.Class <- as.numeric(predict(qda(x = XTrain%*%RP1,  grouping = YTrain), XTest%*%RP1)$class)
    }
    if (base == "knn"){
        if(estmethod == "resub") stop("resubstitution estimate unsuitable for knn classifier") 
        if(estmethod == "loo"){
            for (j in 1:B2) {
                RP <- RPGenerate(p, d, projmethod)
                kcv.voteRP <- sapply(k, function(x) {mean(knn.cv(XTrain%*%RP, YTrain, x) != YTrain, na.rm = TRUE)})
                weight.test <- min(kcv.voteRP)
                if (weight.test <= w) {
                    w <- weight.test
                    RP1 <- RP
                    k1 <- order(kcv.voteRP)[1]
                }
            }
            Train.Class <- as.numeric(knn.cv(XTrain%*%RP1, YTrain, k1))
            Test.Class <- as.numeric(knn(XTrain%*%RP1, XTest%*%RP1, YTrain, k1))
        }
    }
    if (base == "Other") {
        for (j in 1:B2) {
            RP <- RPGenerate(p, d, projmethod)
            weight.test <- mean(Other.classifier(x = XTrain%*%RP, grouping = YTrain, CV = TRUE)$class != YTrain, na.rm = TRUE)
            if (weight.test <= w) {
                w <- weight.test
                RP1 <- RP
            }
        }
        Train.Class <- as.numeric(Other.classifier(x = XTrain%*%RP1, grouping = YTrain, CV = TRUE)$class)
        Test.Class <- as.numeric(Other.classifier(x = XTrain%*%RP1, grouping = YTrain, XTest%*%RP1)$class)
    }
    return(c(Train.Class, Test.Class)) 
}
