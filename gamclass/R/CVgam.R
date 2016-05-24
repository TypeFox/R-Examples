CVgam <-
function (formula, data, nfold = 10, debug.level = 0, method = "GCV.Cp",
              printit = TRUE, cvparts = NULL, gamma = 1, seed = 29)
{
    if (is.null(cvparts)) {
        set.seed(seed)
        cvparts <- sample(1:nfold, nrow(data), replace = TRUE)
    }
    folds <- unique(cvparts)
    khat <- hat <- numeric(nrow(data))
    scale.gam <- summary(gam(formula, data = data, method = method))$scale
    for (i in folds) {
        trainrows <- cvparts != i
        testrows <- cvparts == i
        elev.gam <- gam(formula, data = data[trainrows, ], method = method,
                        gamma = gamma)
        hat[testrows] <- predict(elev.gam, newdata = data[testrows,
                                           ], select = TRUE)
        res <- residuals(elev.gam)
    }
    y <- eval(formula[[2]], envir = as.data.frame(data))
    res <- y - hat
    cvscale <- sum(res^2)/length(res)
    prntvec <- c(GAMscale = scale.gam, `CV-mse-GAM ` = cvscale)
    if (printit)
        print(round(prntvec, 4))
    invisible(list(fitted = hat, resid = res, cvscale = cvscale,
                   scale.gam = scale.gam))
}
