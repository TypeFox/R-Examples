# Copyright 2010 by Roger Bivand and Eric Blankmeyer

lagmess <- function(formula, data = list(), listw, zero.policy=NULL,
    na.action=na.fail, q=10, start=-2.5, control=list(), method="BFGS",
    verbose=NULL) {
    stopifnot(inherits(listw, "listw"))
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
    if (listw$style != "W") warning("weights should be row-stochastic")
    if (class(formula) != "formula") formula <- as.formula(formula)
    mt <- terms(formula, data = data)
    mf <- lm(formula, data, na.action=na.action, method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }
    y <- model.extract(mf, "response")
    stopifnot(all(is.finite(y)))
    X <- model.matrix(mt, mf)
    stopifnot(all(is.finite(X)))
    nullLL <- logLik(lm(formula, data, na.action=na.action))

    W <- as(listw, "CsparseMatrix")
    Y <- powerWeightsMESS(W, y, q=q)

    v <- 0:(q-1)
    G1 <- diag(1/factorial(v))

    bestmess <- optim(start, mymess, gr=NULL, method=method, hessian=TRUE,
        control=control, Y=Y, X=X, G1=G1, v=v)
    alpha <- bestmess$par[1]
    alphase <- 1.0/(bestmess$hessian[1,1])^0.5
    rho <- 1.0 - exp(alpha[1])

    va <- alpha^v
    Sy <- Y %*% G1 %*% va
    data$Sy <- Sy
    formula[[2]] <- formula(~ Sy)[[2]]
    lmobj <- lm(formula=formula, data=data)

    call <- match.call()
    lmobj$call <- call

    res <- list(lmobj=lmobj, alpha=alpha, alphase=alphase, rho=rho,
        bestmess=bestmess, q=q, start=start, na.action=na.act,
        nullLL=nullLL)
    class(res) <- "lagmess"
    res
}

powerWeightsMESS <- function(W, y, q=10) {
        n <- dim(W)[1]
        res <- matrix(NA, nrow=n, ncol=q)
        res[,1] <- y
        last <- W %*% y
        res[,2] <- c(last[,1])
        for (i in 3:q) {
            last <- W %*% last
            res[,i] <- c(last[,1])
        }
        res
}

mymess <- function(alpha, Y, X, G1, v, verbose=FALSE) {
    va <- alpha^v
    Sy <- Y %*% G1 %*% va
    lmobj <- lm(Sy ~ X - 1)
    res <- -c(logLik(lmobj))
    if (verbose) cat("res:", res, "\n")
    res
}

print.lagmess <- function(x, ...) {
    print(x$lmobj, ...)
    cat("Alpha: ", x$alpha, "\n", sep="")
    invisible(x)
}

print.summary.lagmess <- function(x, digits = max(5, .Options$digits - 3),
    signif.stars = FALSE, ...) {
    cat("Matrix exponential spatial lag model:\n")
    print(x$lmsum, signif.stars=signif.stars, digits=digits)
    cat("Alpha: ", format(signif(x$alpha, digits)), ", standard error: ",
        format(signif(x$alphase, digits)), "\n    z-value: ", 
        format(signif((x$alpha/x$alphase), digits)), ", p-value: ",
        format.pval(2 * (1 - pnorm(abs(x$alpha/x$alphase))), digits),
        "\n", sep="")
    res <- x$LR
    cat("LR test value: ", format(signif(res$statistic, digits)), 
        ", p-value: ", format.pval(res$p.value, digits), "\n", sep="")
    cat("Implied rho:", x$rho, "\n")
    cat("\n")
    invisible(x)
}

summary.lagmess <- function(object, ...) {
    object$lmsum <- summary(object$lmobj, ...)
    object$LR <- LR1.lagmess(object)
    class(object) <- "summary.lagmess"
    object
}

LR1.lagmess <- function(object) {
    LLx <- logLik(object)
    LLy <- object$nullLL
    statistic <- 2*(LLx - LLy)
    attr(statistic, "names") <- "Likelihood ratio"
    parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
    if (parameter < 1) 
	stop("non-positive degrees of freedom: no test possible")
    attr(parameter, "names") <- "df"
    p.value <- 1 - pchisq(abs(statistic), parameter)
    estimate <- c(LLx, LLy)
    attr(estimate, "names") <- c("Log likelihood of MESS fit",
        "Log likelihood of OLS fit")
    method <- "Likelihood Ratio diagnostics for spatial dependence"
    res <- list(statistic=statistic, parameter=parameter,
	p.value=p.value, estimate=estimate, method=method)
    class(res) <- "htest"
    res
}

residuals.lagmess <- function(object, ...) {
    object$lmobj$residuals
}

deviance.lagmess <- function(object, ...) {
    deviance(object$lmobj)
}

coef.lagmess <- function(object, ...) {
    ret <- NULL
    ap <- object$alpha
    names(ap) <- "alpha"
    ret <- c(ret, ap)
    ret <- c(ret, coef(object$lmobj))
    ret
}

fitted.lagmess <- function(object, ...) {
    object$lmobj$fitted.values
}

logLik.lagmess <- function (object, ...) 
{
    LL <- c(logLik(object$lmobj))
    class(LL) <- "logLik"
    N <- length(residuals(object))
    attr(LL, "nall") <- N
    attr(LL, "nobs") <- N
    attr(LL, "df") <- object$lmobj$rank + 2
    LL
}

