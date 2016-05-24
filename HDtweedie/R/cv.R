cv.HDtweedie <- function(x, y, group=NULL, p, weights, lambda = NULL, pred.loss = c("deviance", 
    "mae", "mse"), nfolds = 5, foldid, ...) {
    if (missing(p)) p <- 1.5
    if (missing(pred.loss)) 
        pred.loss <- "default" else pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    if (missing(weights)) weights <- rep(1.0,N)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    tweediegrpnet.object <- HDtweedie(x, y, group, p, weights, lambda = lambda, ...)
    lambda <- tweediegrpnet.object$lambda
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
        outlist[[i]] <- HDtweedie(x = x[!which, , drop = FALSE], y = y_sub, group = group, p = p, weights = weights[!which], lambda = lambda, ...)
    }
    ###What to do depends on the pred.loss
    cvstuff <- cv.tweediegrppath(outlist, lambda, x, y, p, weights, foldid, pred.loss)
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
        cvlo = cvm - cvsd, name = cvname, HDtweedie.fit = tweediegrpnet.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.HDtweedie"
    obj
}

cv.tweediegrppath <- function(outlist, lambda, x, y, p, weights, foldid, pred.loss) {
    typenames <- c(deviance = "Tweedie Deviance", mse = "Mean Square Error", mae = "Mean Absolute Error")
    if (pred.loss == "default") 
        pred.loss <- "deviance"
    if (!match(pred.loss, c("deviance", "mse", "mae"), FALSE)) {
        warning("Only 'deviance', 'mse' and 'mae' available; 'deviance' used")
        pred.loss <- "deviance"
    }
    predmat <- matrix(NA, length(y), length(lambda))
	nfolds <- max(foldid)
	nlams <- double(nfolds)
	for (i in seq(nfolds)) {
		which <- foldid==i
		fitobj <- outlist[[i]]
		preds <- predict(fitobj, x[which,], type = "response")
		nlami <- length(outlist[[i]]$lambda)
		predmat[which, seq(nlami)] <- preds
		nlams[i] <- nlami	
	}	
	N <- length(y)
	cvraw <- switch(pred.loss, "deviance"=devi(y, predmat, p), "mae"=abs(y-predmat), "mse"=(y-predmat)^2)
	cvob <- cvcompute(cvraw, weights, foldid, nlams)
	cvraw <- cvob$cvraw
	weights <- cvob$weights
	N <- cvob$N
	cvm <- apply(cvraw, 2, weighted.mean, w=weights, na.rm=TRUE)
	cvsd <- sqrt(apply(scale(cvraw, cvm, scale=FALSE)^2, 2, weighted.mean, w=weights, na.rm=TRUE)/(N-1))
	list(cvm=cvm, cvsd=cvsd, name=typenames[pred.loss])
}


