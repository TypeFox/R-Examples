setMethod("show", "Lda", function(object){

    if(!is.null(cl <- object@call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }

    digits = max(3, getOption("digits") - 3)
    cat("\nPrior Probabilities of Groups:\n")
    print(object@prior)
    if(length(object@center) > 0)
    {
        cat("\nGroup means:\n")
        print(object@center)
        cat("\nWithin-groups Covariance Matrix:\n")
        print(object@cov)
    }
    cat("\nLinear Coeficients:\n")
    print(object@ldf)
    cat("\nConstants:\n")
    print(object@ldfconst)

#    svd <- x$svd
#    names(svd) <- dimnames(x$scaling)[[2]]
#    if(length(svd) > 1) {
#        cat("\nProportion of trace:\n")
#        print(round(svd^2/sum(svd^2), 4), ...)
#    }
    invisible(object)
})


setMethod("predict", "Lda", function(object, newdata){

    ct <- FALSE
    if(missing(newdata))
    {
        newdata <- object@X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)

    if(length(object@center)>0 & ncol(x) != ncol(object@center) | ncol(x) != ncol(object@ldf))
        stop("wrong number of variables")

    ldf <- object@ldf
    ldfconst <- object@ldfconst
    ret <- .mypredictLda(object@prior, levels(object@grp), ldf, ldfconst, x)
    if(ct)
        ret@ct <- mtxconfusion(object@grp, ret@classification)

    ret
})


.mypredictLda <- function(prior, lev, ldf, ldfconst, x){

    ng <- length(prior)
    nm <- names(prior)
    if(is.null(nm))
        nm <- lev

    xx <- x %*% t(ldf)
    xx <- t(t(xx) + ldfconst)
    posterior <- xx
    for(i in 1:nrow(xx)){
        tmp <- sum(exp(xx[i,]))
        for(j in 1:ncol(xx))
            posterior[i,j] <- exp(xx[i,j])/tmp
    }

    cl <- factor(nm[max.col(xx)], levels = lev)
    new("PredictLda", classification=cl, posterior=posterior, x = xx)
}

setMethod("show", "PredictLda", function(object){

    if(!is.null(object@ct))
    {
        tab <- object@ct
        acctab <- t(apply(tab, 1, function(x) x/sum(x)))
        dimnames(acctab) <- dimnames(tab)
        AER <- 1 - sum(diag(tab))/sum(tab)

        prt <- as.matrix(round(c("Apparent error rate" = AER),4))
        colnames(prt) <- ""
        print(prt)

        cat("\nClassification table", "\n")
        print(tab)
        cat("\nConfusion matrix", "\n")
        print(round(acctab, 3))
    }
    else
        print(object@classification)

##    print(object@posterior)
##    print(object@x)
    invisible(object)
})

setMethod("summary", "Lda", function(object, ...){
        new("SummaryLda", ldaobj=object)
})

setMethod("show", "SummaryLda", function(object){
    show(object@ldaobj)
    invisible(object)
})

.wcov.wt <- function (x, gr, wt = rep(1, nrow(x)))
{
    xcall <- match.call()
    if (is.data.frame(x))
        x <- as.matrix(x)
    else if (!is.matrix(x))
        stop("'x' must be a matrix or a data frame")

    if (!all(is.finite(x)))
        stop("'x' must contain finite values only")
    n <- nrow(x)
    p <- ncol(x)
    lev <- levels(as.factor(gr))
    ng <- length(lev)
    dimn <- dimnames(x)


    if(with.wt <- !missing(wt)) {
        if(length(wt) != n)
            stop("length of 'wt' must equal the number of rows in 'x'")
        if(any(wt < 0) || (s <- sum(wt)) == 0)
            stop("weights must be non-negative and not all zero")
    }

    sqrtwt <- sqrt(wt)
    means <- tapply(sqrtwt*x, list(rep(gr, p), col(x)), sum) /rep(tapply(sqrtwt, gr, sum), p)
    wcov <- crossprod(sqrtwt*(x-means[gr,]))/(sum(sqrtwt)-ng)
    dimnames(wcov) <- list(dimn[[2]], dimn[[2]])
    dimnames(means) <- list(lev, dimn[[2]])
    ans <- list(call=xcall, means=means, wcov=wcov, method="Reweighted")
    class(ans) <- "wcov"
    ans
}

print.wcov <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }

    cat("\nGroup means:\n")
    print(x$means, ...)

    cat("\nWithin-groups Covariance Matrix:\n")
    print(x$wcov, ...)

    invisible(x)
}
