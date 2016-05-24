
FADA = function (faobject, K=10,B=20, nbf.cv = NULL,method = c("glmnet", 
    "sda", "sparseLDA"), sda.method = c("lfdr", "HC"), alpha=0.1, ...) 
{
    nbclass <- length(unique(faobject$groups))
    method <- match.arg(method)
    sda.method <- match.arg(sda.method)
    if (!hasArg(maxnbfactors)) {maxnbfactors <- faobject$maxnbfactors}
    if (!hasArg(nfolds)) { nfolds <- faobject$nfolds}
    if (!hasArg(grouped)) {grouped <- faobject$grouped}
    if (!hasArg(min.err)) {min.err <- faobject$min.err}
    if (!hasArg(EM)) {EM <- faobject$EM}
    if (!hasArg(maxiter)) {maxiter <- faobject$maxiter}
    if (! is.null(faobject$fa.testing)) {
        out <- FADA.tmp(faobject= faobject, method, sda.method, ...)
        proba.train <- out$proba.train
        proba.test <- out$proba.test
        predict.test <- out$predict.test
        selected <- out$selected
        out <- out$mod
        cv.error <- NULL
        cv.error.se <- NULL
    }
    else {
        fadta <- faobject$fa.training
        groups <- faobject$groups
        p <- ncol(fadta)
        if (method == "glmnet") {
            out <- LassoML(list(x = fadta, y = groups),...)
            selected <- out$selected
           out <- out$model
           proba.train <- predict(out,fadta,type="response")
        }
        if (method == "sda") {
            ranking.LDA <- sda::sda.ranking(fadta, groups, 
                verbose = FALSE)
            if (sda.method == "lfdr") {selected <- as.numeric(ranking.LDA[ranking.LDA[, "lfdr"] < 0.8, "idx"])}
            if (sda.method == "HC") { thr <- which.max(ranking.LDA[1:round(alpha * p), "HC"]) ;
                selected <- as.numeric(ranking.LDA[1:thr, "idx"]) }
            out <- sda::sda(fadta[, selected,drop=FALSE], groups, verbose = FALSE)
            proba.train <- sda::predict.sda(out,fadta[,selected,drop=FALSE],verbose=FALSE)$posterior
        }
        if (method == "sparseLDA") {
            Xc <- normalize(fadta)
            Xn <- Xc$Xc
            out <- sparseLDA::sda(Xn, factor(groups), ...)
            selected <- out$varIndex
            proba.train <- sparseLDA::predict.sda(out,Xn)$posterior
        }
       cv.out <- crossval(cv.FADA, faobject$data.train,faobject$groups,nbf.cv=nbf.cv,method = method, sda.method = sda.method, alpha = alpha,K=K,B=B,EM=EM,maxiter=maxiter, maxnbfactors= maxnbfactors, min.err= min.err, ...)
        proba.test <- NULL
        predict.test <- NULL
        cv.error <- cv.out$stat
        cv.error.se <- cv.out$stat.se
    }    
    return(list(method = method, selected = selected, proba.train=proba.train,proba.test = proba.test, 
        predict.test = predict.test, cv.error = cv.error,cv.error.se=cv.error.se, mod=out))
}
