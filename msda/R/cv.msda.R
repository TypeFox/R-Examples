cv.msda <- function(x, y, nfolds = 5, lambda = NULL, lambda.opt = "min", 
    ...) {
    y <- drop(y)
    n <- nrow(x)
    p <- ncol(x)
    K <- length(unique(y))
    prior <- rep(0, K)
    for (i in 1:K) {
        prior[i] <- mean(y == i)
    }
    ### Fit the model once to get dimensions etc of output
    tmp <- msda(x, y, lambda = lambda, ...)
    lambda <- tmp$lambda
    ### Now fit the nfold models and store them
    foldid <- sample(rep(seq(nfolds), length = n))
    if (nfolds < 3) 
        stop("nfolds must be bigger than 3; nfolds=5 recommended")
    if (nfolds > n) 
        stop("The number of folds should be smaller than the sample size.")
    residmat <- matrix(NA, nfolds, length(lambda))
    # good <- matrix(0, nfolds, length(lambda))
    for (i in seq(nfolds)) {
        which <- foldid == i
        # cat("\nFold ", i, " is running.\n")
        fit <- msda(x[!which, , drop = FALSE], y[!which], lambda = lambda, 
            ...)
        preds <- predict(fit, x[which, , drop = FALSE])
        nlami <- length(fit$lambda)
        residmat[i, seq(nlami)] <- colMeans(y[which] != preds)
        # good[i, seq(nlami)] <- 1
    }
	residmat[is.na(residmat)] <- 1
    # rN <- colSums(good)
    cvm <- colMeans(residmat, na.rm = TRUE)
    cvsd <- sqrt(colMeans(scale(residmat, cvm, FALSE)^2, na.rm = TRUE)/(nfolds - 
        1))
    if (lambda.opt == "min") {
        lambda.min <- min(lambda[which(cvm == min(cvm, na.rm = TRUE))])
    } else {
        lambda.min <- max(lambda[which(cvm == min(cvm, na.rm = TRUE))])
    }
    id.min <- which.min(cvm)
    lambda.1se <- max(lambda[cvm < min(cvm, na.rm = TRUE) + cvsd[id.min]], 
        na.rm = TRUE)
    # cvm <- na.omit(cvm)
    # cvsd <- na.omit(cvsd)
    # min_len <- min(length(cvm), length(cvsd))
    obj <- list(lambda = lambda, cvm = cvm, cvsd = cvsd, 
        lambda.min = lambda.min, lambda.1se = lambda.1se, msda.fit = tmp)
    class(obj) <- "cv.msda"
    obj
}