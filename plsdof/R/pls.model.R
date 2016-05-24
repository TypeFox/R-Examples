pls.model<-function (X, y, m = ncol(X), Xtest = NULL, ytest = NULL, compute.DoF = FALSE, 
    compute.jacobian = FALSE, use.kernel = FALSE,method.cor="pearson") 
{
    if (compute.DoF == FALSE) {
        compute.jacobian == FALSE
    }
    n <- nrow(X)
    p <- ncol(X)
    X0 <- X
    DoF.max = min(n - 1, p + 1)
    m <- min(m, DoF.max)
    DoF <- NULL
    mse <- cor<-NULL
    prediction <- NULL
    coefficients <- NULL
    sigmahat <- NULL
    RSS <- NULL
    intercept <- NULL
    if (compute.jacobian == TRUE) {
        if (use.kernel == TRUE) {
            pls.object <- kernel.pls.fit(X, y, m, compute.jacobian = compute.jacobian, 
                DoF.max = DoF.max)
        }
        if (use.kernel == FALSE) {
            pls.object <- linear.pls.fit(X, y, m, compute.jacobian = TRUE, 
                DoF.max = DoF.max)
        }
        coefficients <- pls.object$coefficients
        intercept <- pls.object$intercept
        Yhat <- pls.object$Yhat
        yhat <- pls.object$yhat
        RSS <- pls.object$RSS
        sigmahat <- pls.object$sigmahat
        DoF <- pls.object$DoF
    }
    if (compute.jacobian == FALSE) {
        if (use.kernel == TRUE) {
            pls.object <- kernel.pls.fit(X, y, m, compute.jacobian = FALSE, 
                DoF.max = DoF.max)
            intercept <- pls.object$intercept
            coefficients <- pls.object$coefficients
            Yhat <- pls.object$Yhat
            yhat <- pls.object$yhat
            RSS <- pls.object$RSS
            sigmahat <- pls.object$sigmahat
            DoF <- pls.object$DoF
            if (compute.DoF == TRUE) {
                mean.X <- apply(X, 2, mean)
                sd.X <- apply(X, 2, sd)
                sd.X[sd.X == 0] = 1
                X <- X - rep(1, nrow(X)) %*% t(mean.X)
                X <- X/(rep(1, nrow(X)) %*% t(sd.X))
                K <- X %*% t(X)
                dof.object = pls.dof(pls.object, K = K, y = y, 
                  n = n, m = m, DoF.max = DoF.max - 1)
                DoF = c(0, dof.object$DoF) + 1
                sigmahat = c(sqrt(RSS[1]/(n - 1)), dof.object$sigmahat)
            }
        }
        if (use.kernel == FALSE) {
            pls.object <- linear.pls.fit(X, y, m, compute.jacobian = FALSE, 
                DoF.max = DoF.max)
            sigmahat <- pls.object$sigmahat
            DoF <- pls.object$DoF
            Yhat <- pls.object$Yhat
            yhat <- pls.object$yhat
            RSS <- pls.object$RSS
            if (compute.DoF == TRUE) {
                mean.X <- apply(X, 2, mean)
                sd.X <- apply(X, 2, sd)
                sd.X[sd.X == 0] = 1
                X <- X - rep(1, nrow(X)) %*% t(mean.X)
                X <- X/(rep(1, nrow(X)) %*% t(sd.X))
                K <- X %*% t(X)
                dof.object = pls.dof(pls.object, K = K, y = y, 
                  n = n, m = m, DoF.max = DoF.max - 1)
                DoF = c(0, dof.object$DoF) + 1
                sigmahat = c(sqrt(RSS[1]/(n - 1)), dof.object$sigmahat)
            }
            coefficients <- pls.object$coefficients
            intercept <- pls.object$intercept
        }
    }
    if (is.null(Xtest) == FALSE) {
        prediction = rep(1, nrow(Xtest)) %*% t(intercept) + Xtest %*% 
            coefficients
        if (is.null(ytest) == FALSE) {
            res <- matrix(, nrow(Xtest), m + 1)
            for (l in 1:(m + 1)) {
                res[, l] = ytest - prediction[, l]
                if (l>1){
                cor[l]<-cor(ytest,prediction[,l],method=method.cor)
                }
            }
            mse = apply(res^2, 2, mean)
        }
    }
    if (compute.DoF == FALSE) {
        DoF = 1:(m + 1)
        sigmahat = sqrt(RSS/(n - DoF))
    }
    covariance = pls.object$covariance
    return(list(prediction = prediction, mse = mse, cor=cor,coefficients = coefficients, 
        intercept = intercept, DoF = DoF, RSS = RSS, Yhat = Yhat, 
        sigmahat = sigmahat, yhat = yhat, covariance = covariance))
}
