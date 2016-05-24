cv.lspath <- function(outlist, lambda, x, y, foldid, 
    pred.loss, delta) {
    typenames <- c(misclass = "Misclassification Error", loss = "Margin Based Loss")
    if (pred.loss == "default") 
        pred.loss <- "loss"
    if (!match(pred.loss, c("loss"), FALSE)) {
        warning("Only 'loss' available for least squares regression; 'loss' used")
        pred.loss <- "loss"
    }
    ###Turn y into c(0,1)
    y <- as.double(y)
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
    cvraw <- (y-predmat)^2
    cvob <- cvcompute(cvraw, foldid, nlams)
    cvraw <- cvob$cvraw
    N <- cvob$N
    cvm <- apply(cvraw, 2, mean, na.rm = TRUE)
    cvsd <- sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(N - 
        1))
    list(cvm = cvm, cvsd = cvsd, name = typenames[pred.loss])
} 
