# gp <- function(formula, data, 
#                inputs, cov,
#                estim = TRUE, 
#                multistart = NULL,
#                nCores = 1,
#                .export = NULL,
#                .packages = NULL,
#                ...){
#   if (length(multistart)==0) {
#     try(gp1(formula = formula, data = data, 
#         inputs = inputs, cov = cov, estim = estim, ...))
#   } else {
#     nfit <- ncol(multistart)
#     cl <- makeCluster(nCores)
#     registerDoParallel(cl)
#     if (requireNamespace("foreach", quietly = TRUE)){
#       olist <- foreach(i=1:nfit, .errorhandling = 'remove', .verbose=TRUE,
#                        .export = .export, .packages = .packages) %dopar% (
#         gp1(formula = formula, data = data, 
#             inputs = inputs, cov = cov, estim = estim, 
#             parCovIni = multistart[, i],  
#             ...)
#       )
#     }
#     
#     stopCluster(cl)
#     
#     # get the best result
#     nlist <- length(olist)
#     print(nlist)
#     optValue <- vector(length = nlist)
#     if (nlist==0) stop("Maximum Likelihood error")
#     for (i in 1:nlist){
#       optValue[i] <- olist[[i]]$optValue
#     }
#     bestIndex <- which.max(optValue)
#     return(olist[[bestIndex]])
#   }
# }

gp <- function(formula, data,
               inputs, cov,
               estim = TRUE,
               ...) {
    
    Call <- match.call()
    if (! all(inputs %in% colnames(data))) {
        stop("all elements of 'inputs' must be colnames of 'data'")
    }
    if (as.character(formula[[2]]) %in% inputs) {
        stop("the response can not appear in 'inputs'")
    }
    
    X <- as.matrix(data[ , inputs, drop = FALSE], rownames.force = TRUE)
    mf <- model.frame(formula, data = data)
    tt <- terms(mf)
    y <- model.response(mf, "numeric")
    F <- model.matrix(formula, data = mf)
    n <- length(y)
    p <- ncol(F)
    d <- ncol(X)

    if (estim) {
        if (is(cov, "covMan")) {
            ## this is the default value in mle, hence...
            LCall <- as.list(Call)
            if (!("compGrad" %in% names(LCall))) compGrad <- TRUE
            else compGrad <- eval(LCall[["compGrad"]])
            if (compGrad && !cov@hasGrad) {
                stop("when 'compGrad' is given and is TRUE, 'cov' object ",
                     "must compute the gradient")
            }
        }
        
        fit <- try(mle(object = cov,
              y = y, X = X, F = F,
              ...))

        if (inherits(fit, "try-error")) stop("Maximum Likelihood error")
        optValue <- fit$opt$val
        ## replace 'cov'.
        cov <- fit$cov
        ## varNoise <- fit$varNoise
        fit <- fit$trendRes
        
    } else {
        ## if (is.na(varNoise)) varNoise <- NULL
        fit <- try(gls(object = cov,
                       y = y, X = X, F = F,
                       ...))
        if (inherits(fit, "try-error")) stop("Maximum Likelihood error")
        optValue <- NULL
    }
    
    ##   if (is.null(coefTrend)) {
    ##     auxVar <- gls(object = cov, y = y, X = X, F = F,
    ##                   varNoise = varNoise, checkNames = FALSE)
    ##   } 
    res <- list(call = Call,
                terms = tt,
                inputNames = inputs,
                dim = list(n = n, p = p, d = d),
                y = y, X = X, F = F,
                L = fit$L,
                betaHat = fit$betaHat,
                eStar = fit$eStar,
                FStar = fit$FStar,
                sseStar = fit$sseStar,
                RStar = fit$RStar,
                covariance = cov,
                noise = fit$noise,
                varNoise = fit$varNoise,
                optValue = optValue)
    
    class(res) <- "gp"
    return(res)
    
}

##============================================================================
## The 'predict' method is similar to that of DiceKriging::km but does not
## assume stationarity.
##
## TODO: if one wants to have the prediction without covariance matrix,
## a 'diagOnly' argument could be added to the 'covMat' method.
##
## Note that the result is always a vector of pointwise predictions,
## possibly having some attributes.
##
##============================================================================
predict.gp <- function(object, newdata, type, 
                       seCompute = TRUE, covCompute = FALSE,
                       lightReturn = FALSE, biasCorrect = FALSE,
                       forceInterp = FALSE,
                       ...) {
    
    nms <- object$inputNames                 ## MMM method here
    p <- object$dim$p                        ## method here to extract 'p'
    
    ## TODO: allow p = 0 here using suitable 'if' lines
    if (missing(newdata)) {
        XNew <- object$X
        FNew <- object$F
    } else {
        XNew <- newdata[ , nms, drop = FALSE]
        tt <- delete.response(terms(object))
        mf <- model.frame(tt, data = data.frame(newdata))
        FNew <- model.matrix(tt, data = mf)
    }
    yHatTrend <- FNew %*% object$betaHat
    cNew <- covMat(object = object$covariance, X = object$X, Xnew = XNew, 
                   checkNames = FALSE)
    
    if (forceInterp){
        ## force interpolation with a nugget with variance =
        ## 'varNoise'.
        ##
        ##     cov(Z(x), Z(x')) = k(x,x') + varNugget * delta(x, x')
        ## 
        ## so one must add 'varNoise' for each pair (x, x') such that
        ## x = x'
        XMat <- as.matrix(object$X)
        XNewMat <- as.matrix(XNew)
        n1 <- nrow(XMat)
        n2 <- nrow(XNewMat)
        out <- .C("C_covMat1Mat2_WhiteNoise", 
                  as.double(XMat), as.integer(n1),
                  as.double(XNewMat), as.integer(n2), 
                  as.integer(ncol(XMat)),
                  as.double(object$varNoise),
                  ans = double(n1 * n2), PACKAGE = "kergp")
        cWN <- matrix(out$ans, n1, n2)
        cNew <- cNew + cWN
    }
    
    cNewStar <- forwardsolve(object$L, cNew)  
    eNew <- t(cNewStar) %*% object$eStar
    yHat <- yHatTrend + eNew
    
    outputList <- list()
    outputList$trend <- yHatTrend
    outputList$mean <- yHat
    
    if (!lightReturn){
        outputList$c <- cNew
        outputList$cStar <- cNewStar 
    }
    
    if (seCompute) {
        m <- nrow(XNew)
        sd2Total <- rep(NA, m)
        for (i in 1:m){   ## MMM prevoir de l'ecrire en C pour covMan
            sd2Total[i] <- covMat(object = object$covariance, 
                                  X = XNew[i, , drop = FALSE], checkNames = FALSE)
        }
        if (forceInterp){
            sd2Total <- sd2Total + object$varNoise
        }
        
        s2Predict1 <- apply(cNewStar, 2, crossprod)
        s2SK <- pmax(sd2Total - s2Predict1, 0)
        
        ## compute c(x)'*C^(-1)*c(x) for x = newdata
        if (type == "SK"){
            s2Predict <- s2SK
            q95 <- qnorm(0.975)
        } else if (type == "UK") {
            FNewError <-  FNew - t(cNewStar) %*% object$FStar
            FNewError <- forwardsolve(t(object$RStar), t(FNewError))
            s2Predict2 <- apply(FNewError, 2, crossprod)
            s2Predict <- s2SK + s2Predict2
            n <- object$dim$n
            if (biasCorrect) {
                s2Predict <- s2Predict * n / (n - p)
            }
            q95 <- qt(0.975, df = n - p)
        }
        
        s2SK <- as.numeric(s2SK)
        s2Predict <- as.numeric(s2Predict)
        lower95 <- yHat - q95 * sqrt(s2Predict)
        upper95 <- yHat + q95 * sqrt(s2Predict)
        outputList$sdSK <- sqrt(s2SK)
        outputList$sd <- sqrt(s2Predict)
        outputList$lower95 <- lower95
        outputList$upper95 <- upper95  
    }
    
    ## covariance matrix of prediction errors
    if (covCompute) {
        CNew <- covMat(object$covariance, X = XNew, checkNames = FALSE)     ## MMM method here
        covCond <- CNew - crossprod(cNewStar, cNewStar)
        ## if (p > 0L) {
        if (type == "UK"){
            FNewError <-  FNew - t(cNewStar) %*% object$FStar
            FNewError <- forwardsolve(t(object$RStar), t(FNewError))
            covCond <- covCond + crossprod(FNewError, FNewError)
            if (biasCorrect){
                n <- object$dim$n
                covCond <- covCond * n / (n - p)
            }
        }

        if (forceInterp){
            ## enssure that 'covCond' is a matrix, else diag
            ## will not extract the diagonal
            covCond <- as.matrix(covCond)
            diag(covCond) <- diag(covCond) + object$varNoise
        }

        outputList$cov <- covCond
    } 
    
    return(outputList)
}


##==================================================
## Leave-One-Out method (similar to DiceKriging::km) 
##==================================================

influence.gp <- function(model, type = "UK", trend.reestim = TRUE, ...) {
  
    case1 <- ((type == "SK") && (!trend.reestim))
    case2 <- ((type == "UK") && (trend.reestim))
    analytic <- (case1 | case2)
    ## If analytic, we can use Dubrule's formulae
    
    X <- as.matrix(model$X)
    y <- as.matrix(model$y)
    L <- model$L
    n <- nrow(X)
    yhat <- sigma2 <- matrix(NA, n, 1) 
    
    beta <- coef(model, "trend") 
    F <- model$F
    
    if (!analytic)	{	
        
        C <- L%*%t(L)            ## should be improved in future versions  
        
        for (i in 1:n) {
            
            F.i <- F[i, ]
            y.but.i <- y[-i, ]
            F.but.i <- F[-i,]
            C.but.i <- C[-i, -i]
            c.but.i <- C[-i, i]
            L.but.i <- t(chol(C.but.i))					
            x <- forwardsolve(L.but.i, y.but.i)
            M <- forwardsolve(L.but.i, F.but.i)
            M <- as.matrix(M)
            
            if (trend.reestim) {    ## reestimation of beta
                l <- lm(x ~ M - 1)
                beta <- as.matrix(l$coef, ncol = 1)
            }
            z <- x - M%*%beta
            Tinv.c <- forwardsolve(L.but.i, c.but.i)      # only a vector in this case
            y.predict.complement <- t(Tinv.c) %*% z
            y.predict.trend <- F.i %*% beta
            
            y.predict <- y.predict.trend + y.predict.complement
            yhat[i] <- y.predict
            
            sigma2.1 <- crossprod(Tinv.c)
            
            total.sd2 <- C[i, i]
            sigma2[i] <- total.sd2 - sigma2.1
            
            if (type == "UK"){
                T.M <- chol(t(M) %*% M)
                sigma2.mat <- forwardsolve(t(T.M), t(F.i - t(Tinv.c) %*% M))
                sigma2.2 <- apply(sigma2.mat, 2, crossprod)
                sigma2[i] <- sigma2[i] + sigma2.2
            }
            
            sigma2[i] <- max(sigma2[i], 0)
            sigma2[i] <- as.numeric(sigma2[i])
            
        }	## end 'for' loop
        
    } else {
        ## fast computation
        Cinv <- chol2inv(t(L))                       ## cost : n*n*n/3
        if (trend.reestim && (type == "UK")) {
            M <- model$FStar    ##FStar = inv(L)*F      cost : n*n*p     to recompute
            Cinv.F <- Cinv %*% F                     ## cost : 2*n*n*p
            T.M <- chol(crossprod(M))                ## cost : p*p*p/3,  neglected
            aux <- forwardsolve(t(T.M), t(Cinv.F))   ## cost : p*p*n,    neglected
            Q <- Cinv - crossprod(aux)               ## cost : 2*n*n*(p-1/2)
            Q.y <- Q%*%y
            ## Remark:   Q <- Cinv - Cinv.F %*% solve(t(M)%*%M) %*% t(Cinv.F)
            ## direct (not so bad actually)
        } else if ((!trend.reestim) & (type == "SK")){
            Q <- Cinv
            Q.y <- Q %*% (y - F %*% beta)
        } else {      
            stop("This case is not implemented yet")
        }
        sigma2 <- 1 / diag(Q)
        epsilon <- sigma2 * (Q.y)  # cost : n, neglected 
        yhat <- as.vector(y - epsilon)
    }
    
    return(list(mean = as.numeric(yhat),
                sd = as.numeric(sqrt(sigma2))))
    
}


##=========================================
## Plot method (similar to DiceKriging::km) 
##=========================================

plot.gp <- function(x, y, kriging.type = "UK", trend.reestim = TRUE, which = 1:3, ...) {
    
    model <- x
    pred <- influence(model, type=kriging.type, trend.reestim=trend.reestim)
    y <- as.matrix(model$y)
    yhat <- pred$mean
    sigma <- pred$sd
    
    resid <- (y - yhat) / sigma
                                        #par(ask=TRUE)
    xmin <- min(min(yhat), min(y))
    xmax <- max(max(yhat), max(y))
    
    if (!all(is.element(which, c(1, 2, 3)))) {
        warning('Incorrect values of which, default is used instead.')
        which <- 1:3
    }
    
    nFig <- length(which)
    opar <- par(mfrow = c(nFig, 1))
    if (is.element(1, which)) {
        plot(x = y[ , 1], y = yhat,
             xlim = c(xmin, xmax), ylim = c(xmin, xmax),
             xlab = "Exact values", ylab = "Fitted values",
             main = "Leave-one-out", ...)
        lines(x = c(xmin, xmax), y = c(xmin, xmax))
    }
    if (is.element(2, which)) {
        plot(resid, xlab = "Index", ylab = "Standardized residuals",
             main = "Standardized residuals", ...)
    }
    if (is.element(3, which)) {
        qqnorm(resid, main = "Normal QQ-plot of standardized residuals") 
        qqline(resid)
    }
    par(opar)
    
    invisible(pred)

}


##============================================================================
## summary and Co
##============================================================================

summary.gp <- function (object, ...) {
    ans <- object
    ## add some slots or information to the crurrent object before returning
    ## it with the proper S3 class.
    class(ans) <- "summary.gp"
    ans
}

print.summary.gp <-
    function(x, digits = max(3L, getOption("digits") - 3L),
             signif.stars = getOption("show.signif.stars"), ...) {
        
        class(x) <- "gp"
        
        ## call
        cat("\nCall:\n", 
            paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
        
        ## n
        cat("\nNumber of observations:", x$dim$n, "\n")
        
        ## Trend
        cat("\nTrend coef.:\n")
        matTrend <- as.matrix(coef(x, type = "trend"), ncol = 1)
        colnames(matTrend) <- "Value"
        print(matTrend)
        
        ## covariance
        cat(sprintf("\nCovariance whith class \"%s\"\n", class(x$covariance)))
        show(x$covariance)
        cat("\n")
        
        ## noise
        if (x$noise) {
            cat(sprintf("Noise variance: %5.3f\n", x$varNoise))
        }
        
        invisible(x)
        
    }


coef.gp <- function(object, type = c("all", "trend", "covariance", "varNoise"), ...){
    type <- match.arg(type)
    switch(type,
           trend = object$betaHat,
           covariance = coef(object$covariance),
           varNoise = c(varNoise = object$varNoise),
           all = c(object$betaHat, coef(object$covariance), varNoise = object$varNoise))
}

print.gp <- function(x, ...){
    print(x$call)
    cat("\nNumber of observations:", x$dim$n, "\n")
    cat("\nTrend coef.:\n")
    matTrend <- as.matrix(coef(x, type = "trend"), ncol=1)
    colnames(matTrend) <- "Value"
    print(matTrend)
    cat("\n")
    show(x$covariance)
    cat("\nNoise:", coef(x, type = "varNoise"), "\n")  
}


## residuals.gp <-
##   function(object,
##            type = c("working","response", "deviance","pearson", "partial"),
##            ...) {
  
## }
