iC10 <-
function(x, seed=25435) {
    set.seed(seed)
    CN <- x$CN
    Exp <- x$Exp
    if (is.null(CN) & is.null(Exp)) {
        stop("Need to provide at least one source of data\n")
    }
    if (!is.null(CN)) {
        if(!is.null(Exp)) {
            cat("running classifier with copy number and expression...\n")
            Exp <- as.matrix(Exp)
            rownames(x$train.CN) <- paste(rownames(x$train.CN), "CN", sep="_")
            rownames(x$train.Exp) <- paste(rownames(x$train.Exp), "Exp", sep="_")
            rownames(CN) <- paste(rownames(CN), "CN", sep="_")
            rownames(Exp) <- paste(rownames(Exp), "Exp", sep="_")
            data.pamr <- list(x=as.matrix(rbind(x$train.CN,x$train.Exp)), y=x$train.iC10,
                              genenames=c(rownames(x$train.CN), rownames(x$train.Exp)))
            model.train <- pamr.train(data.pamr)
            New <- rbind(CN, Exp)
            tmp <- list(x=New)
            New <- pamr.knnimpute(tmp)$x
            cv <- pamr.cv(model.train, data=data.pamr)
            thr <- cv$threshold[which.min(cv$error)]
            Pred <- pamr.predict(model.train, newx=New,
                                 threshold=thr, type="class")
            Prob <- pamr.predict(model.train, newx=New,
                                 threshold=thr, type="posterior")
            Centroids <- pamr.predict(model.train, newx=New,
                                 threshold=thr, type="cent")
            names(Pred) <- colnames(New)
            cl.type <- "CN+Exp"
        } else {
            cat("running classifier with only copy number...\n")
            if (sum(is.na(CN))>0) {
                ids <- which(apply(CN, 1, function(x) any(is.na(x))))
                x$train.CN <- x$train.CN[-ids,]
                CN <- CN[-ids,]
            }
            rownames(x$train.CN) <- paste(rownames(x$train.CN), "CN", sep="_")
            rownames(CN) <- paste(rownames(CN), "CN", sep="_")
            data.pamr <- list(x=as.matrix(x$train.CN), y=x$train.iC10,
                              genenames=rownames(x$train.CN))
            model.train <- pamr.train(data.pamr)
            New <- CN
            cv <- pamr.cv(model.train, data=data.pamr)
            thr <- cv$threshold[which.min(cv$error)]
            Pred <- pamr.predict(model.train, newx=New,
                                 threshold=thr, type="class")
            Prob <- pamr.predict(model.train, newx=New,
                                 threshold=thr, type="posterior")
            Centroids <- pamr.predict(model.train, newx=New,
                                 threshold=thr, type="cent")
            names(Pred) <- colnames(New)
            cl.type <- "CN"
        }
    } else {
        cat("running classifier with only expression...\n")
        Exp <- as.matrix(Exp)
        New <- as.matrix(Exp)
        rownames(x$train.Exp) <- paste(rownames(x$train.Exp), "Exp", sep="_")
        data.pamr <- list(x=as.matrix(x$train.Exp), y=x$train.iC10, genenames=rownames(x$train.Exp))
        data.pamr <- pamr.knnimpute(data.pamr)
        tmp <- list(x=New)
        New <- pamr.knnimpute(tmp)$x
        ## gene.subset.1 <- apply(x$train.Exp, 1, function(x) !any(is.na(x)))
        ## gene.subset.2 <- apply(New, 1, function(x) !any(is.na(x)))
        ## gene.subset <- which(gene.subset.1  & gene.subset.2)
        ## ## Should we impute?
        model.train <- pamr.train(data.pamr)
        cv <- pamr.cv(model.train, data=data.pamr)
        thr <- cv$threshold[which.min(cv$error)]
        Pred <- pamr.predict(model.train, newx=New,
                             threshold=thr, type="class")
        Prob <- pamr.predict(model.train, newx=New,
                             threshold=thr, type="posterior")
        Centroids <- pamr.predict(model.train, newx=New,
                                  threshold=thr, type="cent")
        names(Pred) <- colnames(New)
        cl.type <- "Exp"
    }
    res <- list(class=Pred, posterior=Prob, centroids=Centroids, fitted=New, map.cn=x$map.cn, 
    map.exp=x$map.exp)
    attr(res, "classifier.type") <- cl.type
    attr(res, "CN.by.feat") <- attr(x, "CN.by.feat")
    attr(res, "Exp.by.feat") <- attr(x, "Exp.by.feat")
    attr(res, "ref") <- attr(x, "ref")
    class(res) <- "iC10"
    res
}
