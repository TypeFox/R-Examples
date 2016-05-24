cv.cpernet <- function(x, y, w = 1.0, lambda = NULL, pred.loss = "loss",
                       nfolds = 5, foldid, tau = 0.8, ...) {
    pred.loss <- match.arg(pred.loss)
    N <- nrow(x)
    y <- drop(y)
    # Fit the model once to get lambda
    cpernet.object <- cpernet(x, y, w, lambda = lambda, tau = tau, ...)
    lambda <- cpernet.object$lambda
    # Obtain active set size
    nz <- coef(cpernet.object, type = "nonzero")
    nz[[1]] <- sapply(nz[[1]], length)
    nz[[2]] <- sapply(nz[[2]], length)
    if (missing(foldid)) {
      foldid <- sample(rep(seq(nfolds), length = N))  
    } else nfolds <- max(foldid)
    if (nfolds < 3) 
      stop("nfolds must be at least 3; nfolds=10 recommended")
    outlist <- vector("list", length = nfolds)
    # Fit the nfolds models
    for (i in seq(nfolds)) {
      whichfold <- (foldid == i)
      y_sub <- y[!whichfold]
      outlist[[i]] <- cpernet(x = x[!whichfold, , drop = FALSE], y = y_sub, 
                              w = w, lambda = lambda, tau = tau, ...)
    }
    # Calculate pred.loss and the model fit
    fun <- paste("cv", class(cpernet.object)[[2]], sep = ".")
    cvstuff <- do.call(fun, list(outlist, lambda, x, y, foldid, pred.loss, w, tau))
    cvm <- cvstuff$cvm
    cvsd <- cvstuff$cvsd
    cvname <- cvstuff$name
    out <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, cvupper = cvm + cvsd, 
                cvlower = cvm - cvsd, nzero = as.list(nz), name = cvname,
                cpernet.fit = cpernet.object)
    lamin <- getmin(lambda, cvm, cvsd)
    obj <- c(out, as.list(lamin))
    class(obj) <- "cv.cpernet"
    obj
} 
