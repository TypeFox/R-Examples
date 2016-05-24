cv.gglasso <- function(x, y, group, lambda = NULL, pred.loss = c("misclass", 
    "loss", "L1", "L2"), nfolds = 5, foldid, delta, ...) {
    if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    if (missing(delta)) 
        delta <- 1
    if (delta < 0) 
        stop("delta must be non-negtive")
    gglasso.object <- gglasso(x, y, group, lambda = lambda, delta = delta, ...)
    lambda <- gglasso.object$lambda
    # predict -> coef
    if (missing(foldid)) 
        foldid <- sample(rep(seq(nfolds), length = N)) else nfolds <- max(foldid)
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=10 recommended")
    outlist <- as.list(seq(nfolds))
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
        which <- foldid == i
        y_sub <- y[!which]
        outlist[[i]] <- gglasso(x = x[!which, , drop = FALSE], y = y_sub, group = group, 
            lambda = lambda, delta = delta, ...)
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(gglasso.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, delta))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
        cvlo = cvm - cvsd, name = cvname, gglasso.fit = gglasso.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.gglasso"
    obj
}



cv.hsvm <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for HHSVM classification; 'misclass' used")
        pred.loss <- "misclass"
    }
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, loss = 2 * hubercls(y * predmat, delta), misclass = (y != 
        ifelse(predmat > 0, 1, -1)))
	N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}


cv.logit <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for logistic regression; 'misclass' used")
        pred.loss <- "misclass"
    }
    prob_min <- 1e-05
    fmax <- log(1/prob_min - 1)
    fmin <- -fmax
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    predmat <- pmin(pmax(predmat, fmin), fmax)
    cvraw <- switch(pred.loss, loss = 2 * log(1 + exp(-y * predmat)), misclass = (y != 
        ifelse(predmat > 0, 1, -1)))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}


cv.sqsvm <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "misclass"
    if (!match(pred.loss, c("misclass", "loss"), FALSE)) {
        warning("Only 'misclass' and 'loss' available for squared SVM classification; 'misclass' used")
        pred.loss <- "misclass"
    }
    ###Turn y into c(0,1)
    y <- as.factor(y)
    y <- c(-1, 1)[as.numeric(y)]
    nfolds <- max(foldid)
    predmat <- matrix(NA, length(y), length(lambda))
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE], type = "link")
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, loss = 2 * ifelse((1 - y * predmat) <= 0, 0, 
        (1 - y * predmat))^2, misclass = (y != ifelse(predmat > 0, 1, -1)))
	N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
}


cv.ls <- function(outlist, lambda, x, y, foldid, pred.loss, delta) {
    typenames <- c(L2 = "Least-Squared loss", L1 = "Absolute loss")
    if (pred.loss == "default") 
        pred.loss <- "L2"
    if (!match(pred.loss, c("L2", "L1"), FALSE)) {
        warning("Only 'L2' and 'L1'  available for least squares models; 'L2' used")
        pred.loss <- "L2"
    }
    predmat <- matrix(NA, length(y), length(lambda))
    nfolds <- max(foldid)
    nlams <- double(nfolds)
    for (i in seq(nfolds)) {
        which <- foldid == i
        fitobj <- outlist[[i]]
        preds <- predict(fitobj, x[which, , drop = FALSE])
        nlami <- length(outlist[[i]]$lambda)
        predmat[which, seq(nlami)] <- preds
        nlams[i] <- nlami
    }
    cvraw <- switch(pred.loss, L2 = (y - predmat)^2, L1 = abs(y - predmat))
    N <- length(y) - apply(is.na(predmat), 2, sum)
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
