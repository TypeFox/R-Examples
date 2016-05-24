`predict.wa` <- function(object, newdata,
                         CV = c("none","LOO","bootstrap", "nfold"),
                         verbose = FALSE, n.boot = 100, nfold = 5,
                         ...) {
    if(missing(newdata))
        return(fitted(object))
    newdata <- as.matrix(newdata)
    if(missing(CV))
        CV <- "none"
    CV <- match.arg(CV)
    ## summaries
    spp.train <- object$n.spp
    spp.fossil <- ncol(newdata)
    n.train <- object$n.samp
    n.fossil <- nrow(newdata)
    n.in.train <- sum(colnames(newdata) %in% names(object$wa.optima))
    Dtype <- object$deshrink
    X <- object$orig.x
    ENV <- object$orig.env
    ## tolerance options from model
    O <- object$options.tol
    useN2 <- object$options.tol$useN2
    ## Doing CV?
    if(identical(CV, "none")) {
        want <- names(object$wa.optima) %in%
        colnames(newdata)
        want <- names(object$wa.optima)[want]
        if(object$tol.dw) {
            pred <- WATpred(newdata[,want], object$wa.optima[want],
                            object$model.tol[want], NROW(newdata[,want]),
                            NCOL(newdata[,want]))
        } else {
            pred <- WApred(newdata[,want], object$wa.optima[want])
        }
        pred <- deshrinkPred(pred, coef(object), type = Dtype)
        ## pred can end up being a 1 col matrix
        ## FIXME: really should check why, but for now, drop
        ## empty dimension
        pred <- drop(pred)
    } else {
        ## CV wanted
        if(identical(CV, "LOO")) {
            loo.pred <- matrix(0, ncol = n.train, nrow = n.fossil)
            ##mod.pred <- length(n.train)
            mod.pred <- numeric(n.train)
            ##useN2 <- object$options.tol$useN2
            want <- names(object$wa.optima) %in% colnames(newdata)
            want <- names(object$wa.optima)[want]
            nr <- NROW(X) - 1
            nr.new <- NROW(newdata)
            nc <- NCOL(X)
            nc.want <- length(want)
            for(i in seq_len(n.train)) {
                if(verbose && ((i %% 10) == 0)) {
                    cat(paste("Leave one out sample", i, "\n"))
                    flush.console()
                }
                wa.optima <- w.avg(X[-i,], ENV[-i])
                ## CV for the training set
                if(object$tol.dw) {
                    tol <- w.tol(X[-i, ], ENV[-i], wa.optima,
                                 useN2 = useN2)
                    ## fix up problematic tolerances
                    tol <- fixUpTol(tol, O$na.tol, O$small.tol,
                                    O$min.tol, O$f, ENV[-i])
                    wa.env <- WATpred(X[-i,], wa.optima, tol,
                                      nr, nc)
                    mod.pred[i] <- WATpred(X[i,,drop=FALSE], wa.optima,
                                           tol, 1, nc)
                } else {
                    wa.env <- WApred(X[-i,], wa.optima)
                    mod.pred[i] <- WApred(X[i,,drop=FALSE], wa.optima)
                }
                deshrink.mod <- deshrink(ENV[-i], wa.env, Dtype)
                wa.env <- deshrink.mod$env
                coefs <- coef(deshrink.mod)
                ## LOO model predictions
                mod.pred[i] <- deshrinkPred(mod.pred[i], coefs,
                                            type = Dtype)
                ## newdata predictions
                pred <- if(object$tol.dw) {
                    WATpred(newdata[,want], wa.optima[want],
                            tol[want], nr.new, nc.want)
                } else {
                    WApred(newdata[,want], wa.optima[want])
                }
                loo.pred[,i] <- deshrinkPred(pred, coefs, type = Dtype)
            }
            ## average the LOO predictions
            pred <- rowMeans(loo.pred)
        } else if(identical(CV, "bootstrap")) {
            boot.pred <- matrix(0, ncol = n.boot, nrow = n.fossil)
            oob.pred <- matrix(NA, ncol = n.boot, nrow = n.train)
            want <- names(object$wa.optima) %in% colnames(newdata)
            want <- names(object$wa.optima)[want]
            nr.new <- NROW(newdata)
            nc <- NCOL(X)
            nc.want <- length(want)
            for(i in seq_len(n.boot)) {
                if(verbose && ((i %% 100) == 0)) {
                    cat(paste("Bootstrap sample", i, "\n"))
                    flush.console()
                }
                ## bootstrap sample
                sel <- sample.int(n.train, n.train, replace = TRUE)
                nr <- NROW(X[sel, , drop = FALSE]) ## number of samples
                nr.oob <- NROW(X[-sel, , drop = FALSE])
                wa.optima <- w.avg(X[sel,,drop = FALSE], ENV[sel])
                ## CV for the training set
                if(object$tol.dw) {
                    tol <- w.tol(X[sel, , drop = FALSE], ENV[sel],
                                 wa.optima, useN2 = useN2)
                    ## fix up problematic tolerances
                    tol <- fixUpTol(tol, O$na.tol, O$small.tol,
                                    O$min.tol, O$f, ENV[sel])
                    wa.env <- WATpred(X[sel, , drop = FALSE],
                                      wa.optima, tol, nr, nc)
                    pred <- WATpred(X[-sel, ,drop=FALSE], wa.optima,
                                    tol, nr.oob, nc)
                } else {
                    wa.env <- WApred(X[sel, ,drop = FALSE], wa.optima)
                    pred <- WApred(X[-sel, ,drop=FALSE], wa.optima)
                }
                deshrink.mod <- deshrink(ENV[sel], wa.env, Dtype)
                wa.env <- deshrink.mod$env
                coefs <- coef(deshrink.mod) #$coef
                ## sample specific errors or model performance stats
                oob.pred[-sel,i] <- deshrinkPred(pred, coefs, type = Dtype)
                ## do the prediction step
                pred <- if(object$tol.dw) {
                    WATpred(newdata[,want], wa.optima[want],
                            tol[want], nr.new, nc.want)
                } else {
                    WApred(newdata[,want], wa.optima[want])
                }
                boot.pred[,i] <- deshrinkPred(pred, coefs, type = Dtype)
            }
            pred <- rowMeans(boot.pred)
        } else if (identical(CV, "nfold")) {
            boot.pred <- matrix(0, ncol = n.boot, nrow = n.fossil)
            oob.pred <- matrix(NA, ncol = n.boot, nrow = n.train)
            want <- names(object$wa.optima) %in% colnames(newdata)
            want <- names(object$wa.optima)[want]
            nr.new <- NROW(newdata)
            nc <- NCOL(X)
            nc.want <- length(want)
            ind <- rep(1:nfold, length = n.train)
            for(i in seq_len(n.boot)) {
                if(verbose && ((i %% 100) == 0)) {
                    cat(paste("n-fold sample", i, "\n"))
                    flush.console()
                }
                ## n-fold sample
                pind <- sample(ind)
                for (k in seq_len(nfold)) {
                    sel <- pind != k
                    nr <- NROW(X[sel, , drop = FALSE]) ## number of samples
                    nr.oob <- NROW(X[-sel, , drop = FALSE])
                    wa.optima <- w.avg(X[sel,], ENV[sel])
                    ## CV for the training set
                    if(object$tol.dw) {
                        tol <- w.tol(X[sel, , drop = FALSE], ENV[sel],
                                     wa.optima, useN2 = useN2)
                        ## fix up problematic tolerances
                        tol <- fixUpTol(tol, O$na.tol, O$small.tol,
                                        O$min.tol, O$f, ENV[sel])
                        wa.env <- WATpred(X[sel, , drop = FALSE],
                                          wa.optima, tol, nr, nc)
                        pred <- WATpred(X[-sel, ,drop=FALSE], wa.optima,
                                        tol, nr.oob, nc)
                    } else {
                        wa.env <- WApred(X[sel, ,drop = FALSE], wa.optima)
                        pred <- WApred(X[-sel, ,drop=FALSE], wa.optima)
                    }
                    ## do the model bits
                    #ones <- rep(1, length = length(wa.optima))
                    #miss <- is.na(wa.optima)
                    #ones[miss] <- 0
                    #wa.optima[miss] <- 0
                    #rowsum <- X[sel,] %*% ones
                    #wa.env <- (X[sel,] %*% wa.optima) / rowsum
                    deshrink.mod <- deshrink(ENV[sel], wa.env, Dtype)
                    wa.env <- deshrink.mod$env
                    coefs <- coef(deshrink.mod) #$coef
                    ## if we want sample specific errors or
                    ## model performance stats
                    #rowsum <- X[!sel,] %*% ones
                    #pred <- (X[!sel,] %*% wa.optima) / rowsum
                    oob.pred[!sel,i] <- deshrinkPred(pred, coefs,
                                                     type = Dtype)
                    ## do the prediction step
                    #want <- names(wa.optima) %in% colnames(newdata)
                    #want <- names(wa.optima)[want]
                    #ones <- rep(1, length = length(want))
                    #miss <- miss[want]
                    #ones[miss] <- 0
                    #rowsum <- newdata[,want] %*% ones
                    #pred <- (newdata[,want] %*% wa.optima[want]) /
                    #    rowsum
                    pred <- if(object$tol.dw) {
                        WATpred(newdata[,want], wa.optima[want],
                                tol[want], nr.new, nc.want)
                    } else {
                        WApred(newdata[,want], wa.optima[want])
                    }
                    boot.pred[,i] <- deshrinkPred(pred, coefs, type = Dtype)
                }
            }
            pred <- rowMeans(boot.pred)
        }
    }
    .call <- match.call()
    .call[[1]] <- as.name("predict")
    names(pred) <- rownames(newdata)
    retval <- list(pred = list(pred = pred, rmsep = NULL),
                   performance = NULL, model.pred = NULL)
    if(identical(CV, "none")) {
        retval$performance <- with(object,
                                   list(r.squared = r.squared,
                                        avg.bias = avg.bias,
                                        max.bias = max.bias,
                                        rmsep = rmse))
        retval$model.pred <- list(pred = fitted(object))
    } else if(identical(CV, "LOO")) {
        mod.r.squared <- cor(mod.pred, ENV)^2
        mod.resid <- ENV - mod.pred ##mod.pred - ENV
        mod.avg.bias <- mean(mod.resid)
        mod.max.bias <- maxBias(mod.resid, ENV)
        mod.rmsep <- sqrt(mean(mod.resid^2))
        retval$performance <- list(r.squared = mod.r.squared,
                                   avg.bias = mod.avg.bias,
                                   max.bias = mod.max.bias,
                                   rmsep = mod.rmsep)
        names(mod.pred) <- names(mod.resid) <- rownames(X)
        retval$model.pred <- list(pred = mod.pred,
                                  resid = mod.resid)
    } else {
        mod.pred <- rowMeans(oob.pred, na.rm = TRUE)
        mod.resid <- ENV - mod.pred ##mod.pred - ENV
        s1 <- apply(oob.pred, 1, sd, na.rm = TRUE)
        ##s2 <- sqrt(rowMeans((oob.pred - ENV)^2, na.rm = TRUE))
        s2 <- sqrt(rowMeans((ENV - oob.pred)^2, na.rm = TRUE))
        mod.s1 <- sqrt(mean(s1^2))
        mod.s2 <- sqrt(mean(mod.resid^2))
        samp.rmsep <- sqrt(s1^2 + mod.s2^2)
        mod.rmsep <- sqrt(mod.s1^2 + mod.s2^2)
        mod.r.squared <- cor(mod.pred, ENV)^2
        mod.avg.bias <- mean(mod.resid)
        mod.max.bias <- maxBias(mod.resid, ENV)
        retval$performance <- list(r.squared = mod.r.squared,
                                   avg.bias = mod.avg.bias,
                                   max.bias = mod.max.bias,
                                   rmsep = mod.rmsep)
        names(mod.pred) <- names(mod.resid) <-
            names(samp.rmsep) <- rownames(X)
        retval$model.pred <- list(pred = mod.pred,
                                  resid = mod.resid,
                                  rmsep = samp.rmsep)
        test.s1 <- apply(boot.pred, 1, sd)
        test.rmsep <- sqrt(test.s1^2 + mod.s2^2)
        names(test.rmsep) <- names(retval$pred$pred)
        retval$pred$pred <- pred
        retval$pred$rmsep <- test.rmsep
    }
    retval$call = .call
    if (identical(CV, "nfold"))
        CV <- paste(nfold, "fold", sep="-")
    retval$CV.method <- CV
    retval$deshrink <- Dtype
    retval$tol.dw <- object$tol.dw
    class(retval) <- "predict.wa"
    retval
}

WApred <- function(X, optima) {
    ones <- rep.int(1, length(optima))
    miss <- is.na(optima)
    ones[miss] <- 0
    optima[miss] <- 0
    rsum <- X %*% ones
    ((X %*% optima) / rsum)
}

WATpred <- function(X, optima, tol, nr, nc) {
    miss <- is.na(optima)
    optima[miss] <- 0
    tol[miss] <- 1
    tol2 <- tol^2
    res <- .C("WATpred", spp = as.double(X), opt = as.double(optima),
              tol2 = as.double(tol2), nr = as.integer(nr),
              nc = as.integer(nc), want = as.integer(miss),
              stat = double(nr),
              NAOK = TRUE, PACKAGE="analogue")$stat
    res
}
