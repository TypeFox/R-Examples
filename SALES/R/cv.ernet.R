cv.ernet <- function(x, y, lambda = NULL, pred.loss = "loss", 
            nfolds = 5, foldid, tau = 0.5, ...) {
    pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    ###Fit the model once to get dimensions etc of output
    y <- drop(y)
    ernet.object <- ernet(x, y, lambda = lambda, tau = tau, ...)
    lambda <- ernet.object$lambda
    # predict -> coef
    nz <- sapply(coef(ernet.object, type = "nonzero"), length)
    if (missing(foldid)) {
      foldid <- sample(rep(seq(nfolds), length = N))  
    } else nfolds <- max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be at least 3; nfolds=10 recommended")
    outlist <- vector("list", length = nfolds)
    ###Now fit the nfold models and store them
    for (i in seq(nfolds)) {
      whichfold <- foldid == i
      y_sub <- y[!whichfold]
      outlist[[i]] <- ernet(x = x[!whichfold, , drop = FALSE], 
          y = y_sub, lambda = lambda, tau = tau, ...)
    }
    ###What to do depends on the pred.loss and the model fit
    fun <- paste("cv", class(ernet.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, 
        pred.loss, tau))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + 
        cvsd, cvlower = cvm - cvsd, nzero = nz, name = cvname, ernet.fit = ernet.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.ernet"
    obj
} 
