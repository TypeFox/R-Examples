splitMatrix <- function (x, f, drop=FALSE) {
    lapply(split(seq_len(nrow(x)), f, drop = drop),
           function(ind) x[ind, , drop = FALSE])
}

fixEvent <- function(event)
{
    ### Convert outcome in clogit model to 0/1 binary coding
    
    if (any(is.na(event)))
        stop("Event contains missing values")

    if (is.logical(event)) {
        status <- is.numeric(event)
    }
    else if (is.numeric(event)) {
        status <- if (max(event) == 2) event - 1  else event
        temp <- (status == 0 | status == 1)
        if (!all(temp)) {
            warning("If outcome is numeric then it must be coded 0/1 or 1/2")
        }
    }
    else if (is.factor(event)) {
        if (nlevels(event) != 2)
            stop("If outcome is a factor then it must have 2 levels")
        status <- event == levels(event)[2]
    }
    return(as.integer(status))
}  

isInformative <- function(Xsplit, ysplit, strata)
{
    ## Identify which observations are informative in a conditional
    ## logistic regression.
    
    is.homogeneous <- function(x) { all(x==x[1]) }
    y.bad <- sapply(ysplit, is.homogeneous)
    X.bad <- sapply(Xsplit, function(x) { all(apply(x, 2, is.homogeneous)) })

    is.informative <- vector("list", length(ysplit))
    for (i in seq(along=is.informative)) {
        canuse <- (!y.bad[i]) && (!X.bad[i])
        is.informative[[i]] <- rep(canuse, length(ysplit[[i]]))
    }
    return(unsplit(is.informative, strata))
}
        
fitClogit <- function(X, y, offset, strata, init, iter.max, eps, toler.chol)
{
    ## Safe wrapper around the C function "clogit" that ensures all
    ## arguments have the correct type and storage mode.
    
    y <- fixEvent(y)
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    if (!is.double(X)) {
        X <- matrix(as.double(X), nrow(X), ncol(X))
    }
    if (is.null(offset)) {
        offset <- rep(0, nrow(X))
    }
    offset <- as.double(offset);
    
    ## Split into strata
    Xsplit <- splitMatrix(X, strata, drop=TRUE)
    ysplit <- split(y, strata, drop=TRUE)
    osplit <- split(offset, strata, drop=TRUE)

    info <- isInformative(Xsplit, ysplit, strata)
    if (!any(info)) {
        stop("There are no informative observations")
    }
    
    ans <- .Call("clogit", Xsplit, ysplit, osplit, as.double(init),
                 as.integer(iter.max), as.double(eps), as.double(toler.chol),
                 PACKAGE="Epi")
    ans$informative <- info
    return(ans)
}

clogistic <- function (formula, strata, data, subset, na.action,
                       init, model = TRUE, x = FALSE, y = TRUE,
                       contrasts = NULL, iter.max=20, eps=1e-6,
                       toler.chol = sqrt(.Machine$double.eps)) 
{
    ## User interface, edited version of glm
    
    call <- match.call()
    if (missing(data)) 
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset",
                 "na.action", "offset", "strata"), names(mf), 0L)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (is.null(Y)) stop("missing outcome")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm)) 
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt)) 
        model.matrix(mt, mf, contrasts)
    else stop("Invalid model matrix")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y)) 
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                length(offset), NROW(Y)), domain = NA)
    }

    strata <- model.extract(mf, "strata")
    if (is.null(strata)) stop("argument 'strata' missing")
    
    contrasts <- attr(X, "contrasts")
    if (attr(mt, "intercept") > 0) {
        X <- X[,-1, drop=FALSE]
    }
    if (missing(init))
        init <- rep(0, ncol(X))

    if (iter.max < 0)
        stop("'iter.max' must be non-negative")
    if (eps <= 0)
        stop("'eps' must be positive")
    if (toler.chol <= 0)
        stop("'toler.chol' must be positive")
    if (eps < toler.chol)
        stop("'toler.chol' must be smaller than 'eps'")
    
    fit <- fitClogit(X = X, y = Y, offset = offset, strata=strata, init=init,
                     toler.chol=toler.chol, eps=eps, iter.max=iter.max)
    if (fit$flag <= 0) {
        stop("Information matrix is not positive definite")
    }
    else if (fit$flag == 1000) {
        warning("Iteration limit exceeded")
    }

    nvar <- length(init)
    which.sing <- if (fit$flag < nvar) {
        diag(fit$var)==0
    } else {
        rep(FALSE, nvar)
    }
    fit$coefficients[which.sing] <- NA
    fit$flag <- NULL
    
    ## Add back in parameter names
    cfnames <- colnames(X)
    names(fit$coefficients) <- cfnames
    dimnames(fit$var) <- list(cfnames, cfnames)

    fit$n <- sum(fit$informative)
    if (model) {
        fit$model <- mf
    }
    else {
        ## Without model frame this cannot be interpreted
        fit$informative <- NULL
    }
    fit$na.action <- attr(mf, "na.action")
    if (x) 
        fit$x <- X
    if (!y) 
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt, 
                       contrasts = contrasts, xlevels = .getXlevels(mt, mf)))
    class(fit) <- c("clogistic")
    fit
}

coef.clogistic <- function(object,...) { object$coefficients }

vcov.clogistic <- function(object, ...) { object$var }

print.clogistic <- function (x, digits = max(options()$digits - 4, 3), ...) 
{
    ## Print method for clogistic objects, edited from print.coxph
    
    cat("\nCall: ", deparse(x$call), "\n\n", sep="\n")

    savedig <- options(digits = digits)
    on.exit(options(savedig))
    coef <- coef.clogistic(x)
    se <- sqrt(diag(vcov.clogistic(x)))
    if (is.null(coef) | is.null(se)) 
        stop("Input is not valid")

    coefmat <- cbind(coef, exp(coef), se, coef/se,
                     signif(1 - pchisq((coef/se)^2, 1), digits - 1))
    dimnames(coefmat) <- list(names(coef),
                              c("coef", "exp(coef)", "se(coef)", "z", "p"))
    cat("\n")
    prmatrix(coefmat)
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    if (is.null(x$df)) 
        df <- sum(!is.na(coef))
    else df <- round(sum(x$df), 2)
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), 
        "  on ", df, " df,", " p=", format(1 - pchisq(logtest, df)),
        ", n=", x$n, sep = "")
    cat("\n")
    invisible()
}
