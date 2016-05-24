bsnOpt <-
function (X= matrix(rnorm(25*10), ncol=10), y=NULL,
              method = "exhaustive", nvmax = NULL,
              nbest=1, intercept=TRUE, criterion="cp", tcrit=NULL,
              print.summary = TRUE, really.big = FALSE, ...)
{
    leaps.out <- try(requireNamespace("leaps"), silent = TRUE)
    if ((is.logical(leaps.out) == TRUE) & (leaps.out == TRUE)) {
        if(is.data.frame(X)){
            if(intercept) X <- model.matrix(~., data=X)[,-1] else
            X <- model.matrix(~-1+., data=X)
        }
        m <- dim(X)[1]
        n <- dim(X)[2]
        if (is.null(colnames(X)))
            colnames(X) <- paste("V", 1:ncol(X), sep = "")
        if(is.null(y))y <- rnorm(m)
        u <- leaps::regsubsets(X, y, method = method, nvmax = nvmax,
                        nbest = nbest, really.big = really.big, intercept=intercept, ...)
        usum <- summary(u)
        nmodels <- nrow(usum$which)
        crit <- usum[[criterion]]
        b <- coef(u, id=1:nmodels, vcov=TRUE)
        if(intercept)
            abstmin <- sapply(b,function(x)min(abs(x[-1])/sqrt(diag(attr(x,"vcov"))[-1])))
        else
            abstmin <- sapply(b,function(x)min(abs(x)/sqrt(diag(matrix(attr(x,"vcov"))))))
        numx <- sapply(b, length) - as.numeric(intercept)
        if(!is.null(tcrit)){
            tcheck <- sapply(b,function(x)abstmin>tcrit)
        } else tcheck <- rep(TRUE, nmodels)
        if(any(tcheck)){
            if(criterion%in%c("bic","cp"))
                keep <- which(tcheck)[which.min(crit[tcheck])] else
            if(criterion=="adjr2")
                keep <- which(tcheck)[which.max(crit[tcheck])] else
            stop(paste("Criterion", criterion, "not available"))
            choosecols <- summary(u)$which[keep,]
            colnam <- names(choosecols)[choosecols]
            if(intercept){
                x <- X[, choosecols[-1], drop=FALSE]
                colnames(x) <- colnam[-1]
                u1 <- lm(y ~ x)
            } else {
                x <- X[, choosecols, drop=FALSE]
                colnames(x) <- colnam
                u1 <- lm(y ~ -1+x)
            }
            if (print.summary){
                print(summary(u1, corr = FALSE))
                if(length(colnames(x))==1)
                    cat(paste("NB: x-variable is", " '", colnames(x), "' ", sep=""), "\n")
            }
            invisible(list(best=u1, abstmin=abstmin, regsubsets_obj=u))
        } else invisible(list(best=NULL, abstmin=abstmin, regsubsets_obj=u))
    }
    else {
        print("Error: package leaps is not installed properly")
    }
}
