summary.hyperblm <- function(object, hessian = FALSE,
                             nboots = 1000, ...) {

    if (! "hyperblm" %in% class(object))
        stop("Object must be of class hyperblm")
    x <- object$xMatrix
    y <- as.numeric(object$yVec)
    distparam <- as.numeric(object$distributionParams)
    coef <- as.numeric(object$coef)
    param <- c(distparam, coef)
    hs <- NULL

    if(hessian == FALSE){
        n <- ncol(x)
        bCoef <- matrix(ncol = n, nrow = nboots)

        for(i in 1:nboots){
            w <- apply(x, 2, sample, replace = TRUE)
            z <- as.vector(w%*%as.matrix(coef)) +
                rhyperb(length(y), param = distparam)
            tryOpt <- try(fittedresult <-
                          hyperblmFit(x = w, y = z, method = object$method,
                                      startMethod = object$startMethod,
                                      startStarts = object$startStarts,
                                      paramStart = object$paramStart, ...),
                          silent = TRUE)
            if (class(tryOpt) == "try-error"){
                bCoef[i, ] <- rep(NA, n)
            } else {
                bCoef[i, ] <- as.numeric(fittedresult$coef)
            }
        }
        ses <- sqrt(apply(bCoef, 2, var, na.rm = TRUE))
    } else {
        param <- c(distparam, coef)
        llfunc <- function(param){
            resids <- y - as.vector(x %*% as.matrix(param[-(1:4)]))
            hyperbDens <- dhyperb(resids, param = param[1:4])
            ##cat("log-likelihood is", -sum(log(hyperbDens)), "\n")
            return(sum(log(hyperbDens)))
        }
        hs <- tsHessian(param, llfunc)
        varcov <- solve(hs)
        variance <- diag(varcov)
        variance <- ifelse(abs(variance) < 1e-03, abs(variance), variance)
        ses <- sqrt(variance)[-(1:4)]
    }

    object$tval <- coef/ses
    object$rdf <- nrow(object$xMatrix) - ncol(object$xMatrix) -3
    object$pval <- 2*pt(abs(object$tval), object$rdf, lower.tail = FALSE)
    object$hessian <- hs
    object$sds <- ses

    class(object) <- "summary.hyperblm"
    return(object)
} ## End of summary.hyperblm

### Print summary
print.summary.hyperblm <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{

    if (class(x) != "summary.hyperblm")
        stop("Object must belong to class summary.hyperblm")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    cat("\nData:     ", x$xNames, "\n")


    if (!is.null(x$hessian)) {
        cat ("Hessian: (mu, delta, pi, zeta) parameter\n")
        print.default(x$hessian)
    }
    cat("Parameter estimates:\n")
    if (is.null(x$sds)) {
        coefficients <- x$coef
        names(coefficients) <- x$xNames
        print.default(format(coefficients, digits = digits),
                      print.gap = 2, quote = FALSE, right = TRUE)
    } else {
        coefficients <- cbind(x$coef, x$sds,
                              x$tval, x$pval)
        dimnames(coefficients) <- list(x$xNames,
                                       c("Estimate", "Std. Error",
                                         "t value", "Pr(>|t|)"))
        printCoefmat(coefficients)
    }
    distributionParams <- x$distributionParams
    names(distributionParams) <- c("mu", "delta", "alpha", "beta")
    cat("\nParameters for hyperbolic error distribution:\n")
    print.default(format(distributionParams, digits = digits),
                  print.gap = 2, quote = FALSE, right = TRUE)
    cat("Likelihood:        ", x$mle, "\n")
    cat("Method:            ", x$method, "\n")
    cat("Convergence code:  ", x$conv, "\n")
    cat("Iterations:        ", x$iter, "\n")
    invisible(x)
}




