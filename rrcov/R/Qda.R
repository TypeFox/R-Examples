setMethod("show", "Qda", function(object){

    if(!is.null(cl <- object@call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }

    digits = max(3, getOption("digits") - 3)
    cat("\nPrior Probabilities of Groups:\n")
    print(object@prior)
    cat("\nGroup means:\n")
    print(object@center)

    ng <- length(object@prior)
    for(i in 1:ng){
        cat("\nGroup: ",levels(object@grp)[i], "\n")
        print(object@cov[,,i])
    }

#    cat("\nLinear Coeficients:\n")
#    print(object@ldf)
#    cat("\nConstants:\n")
#    print(object@ldfconst)

#    svd <- x$svd
#    names(svd) <- dimnames(x$scaling)[[2]]
#    if(length(svd) > 1) {
#        cat("\nProportion of trace:\n")
#        print(round(svd^2/sum(svd^2), 4), ...)
#    }
    invisible(object)
})


setMethod("predict", "Qda", function(object, newdata){

    ct <- FALSE
    if(missing(newdata))
    {
        newdata <- object@X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }

    x <- as.matrix(newdata)

    if(ncol(x) != ncol(object@center))
        stop("wrong number of variables")


    ret <- .mypredictQda(object@prior, levels(object@grp), object@center, object@covinv, object@covdet, x)
    if(ct)
        ret@ct <- mtxconfusion(object@grp, ret@classification)

    ret
})


.mypredictQda <- function(prior, lev, center, covinv, covdet, x){

    ng <- length(prior)
    nm <- names(prior)
    if(is.null(nm))
        nm <- lev

    xx <- matrix(0, nrow=nrow(x), ncol=ng)
    posterior <- xx
    for(j in 1:nrow(x)){
        for(i in 1:ng){
            xx[j,i] <- (x[j,]-center[i,]) %*% covinv[,,i] %*% (x[j,]-center[i,]) + log(covdet[i]) - 2*log(prior[i])
            xx[j,i] <- -0.5*xx[j,i]
        }
    }

    for(i in 1:nrow(xx)){
        tmp <- sum(exp(xx[i,]))
        for(j in 1:ncol(xx))
            posterior[i,j] <- exp(xx[i,j])/tmp
    }

    cl <- factor(nm[max.col(xx)], levels = lev)
    new("PredictQda", classification=cl, posterior=posterior, x = xx)
}

setMethod("show", "PredictQda", function(object){

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

setMethod("summary", "Qda", function(object, ...){
        new("SummaryQda", qdaobj=object)
})

setMethod("show", "SummaryQda", function(object){
    show(object@qdaobj)
    invisible(object)
})
