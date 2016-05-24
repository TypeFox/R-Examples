setMethod("show", "Simca", function(object){

    if(!is.null(cl <- object@call)) {
        names(cl)[2] <- ""
        cat("Call:\n")
        dput(cl)
    }

    digits = max(3, getOption("digits") - 3)
    cat("\nPrior Probabilities of Groups:\n")
    print(object@prior)

    cat("\nPca objects for Groups:\n")

    for(i in 1:length(object@pcaobj))
    {
        print(summary(object@pcaobj[[i]]))
    }
    invisible(object)
})


.testmodel <- function(object, newdata)
{

    if(missing(newdata) || is.null(newdata))
    {
        ## --- No new data provided
        stop("'newdata' not provided - cannot continue")
        newX <- model.matrix(object)
    }
    else if(is.matrix(newdata) | is.data.frame(newdata))
    {
        if(is.data.frame(newdata))
            newdata <- data.matrix(newdata)

        ## --- 'newdata' is a matrix - check for proper number of variables
        if(ncol(newdata) != length(getCenter(object@pcaobj[[1]])))
            stop("'newdata' does not have the correct number of columns")
        newX <- newdata
    }
    else
    {
        ## --- 'newdata' is a data.frame
        stop("'newdata' is a model.frame - currently not supported!")
##        Terms <- delete.response(terms(object))
##        m <- model.frame(Terms, newdata, na.action = na.action)
##        if(!is.null(cl <- attr(Terms, "dataClasses")))
##            .checkMFClasses(cl, m)
##        newX <- delete.intercept(model.matrix(Terms, m))
    }

    ng <- length(object@pcaobj)
    n <- nrow(newX)

    sd <- sdsc <- matrix(NA, nrow=n, ncol=ng)
    od <- odsc <- matrix(NA, nrow=n, ncol=ng)
    out <- object@pcaobj
    for(jg in 1:ng)
    {
        for(index in 1:n)
        {
            out[jg] <- object@pcaobj[jg]                              # a Pca object
            dataicentered <- newX[index,] - getCenter(out[[jg]])      # center the data
            scorei <- dataicentered %*% getLoadings(out[[jg]])        # and compute the PCA scores
            dataitilde <- scorei %*% t(getLoadings(out[[jg]]))        # ??

            ev <- getEigenvalues(out[[jg]])
            sd[index, jg] <- sqrt(scorei %*% (diag(1/ev, nrow=length(ev))) %*% t(scorei))
            od[index,jg] <- rrcov::vecnorm(dataicentered - dataitilde)
            if(out[[jg]]@cutoff.od != 0)
            {
                odsc[index, jg] <- od[index, jg] / out[[jg]]@cutoff.od
            }else
            {
                odsc[index, jg] <- 0
            }

            sdsc[index, jg] <- sd[index, jg] / out[[jg]]@cutoff.sd
        }
    }

    ## return the orthogonal and score distances
    list(odsc=odsc, sdsc=sdsc)
}

##  Obtain the group assignments for given od's and sd's.
.assigngroup <- function(odsc, sdsc, method=1, gamma=0.5)
{
    if(method == 1)
    {
        sd <- sdsc
        od <- odsc
    } else if(method == 2)
    {
        sd <- sdsc ^ 2
        od <- odsc ^ 2
    }

    result <- matrix(NA, nrow=nrow(od), ncol=length(gamma))
    for(igamma in 1:length(gamma))
    {
        tdist <- gamma[igamma] * od + (1-gamma[igamma]) * sd     # gamma(igamma) .* od + (1-gamma(igamma)) .* sd
                                                                 # .* in MATLAB is array multiplication. A.*B is the
                                                                 #  element-by-element product of the arrays A and B.
                                                                 #  A and B  must have the same size, unless one of
                                                                 #  them is a scalar.
        for(i in 1:nrow(od))
        {
            result[i, igamma] <- which(tdist[i,] == min(tdist[i,]))
        }
    }

    result
}

setMethod("predict", "Simca", function(object, newdata, method=2, gamma=0.5, ...){

    ct <- FALSE
    if(missing(newdata))
    {
        newdata <- object@X         # use the training sample
        ct <- TRUE                  # perform cross-validation
    }else if(length(dim(newdata)) != 2L)
        stop("'newdata' must be a matrix or data frame")

    x <- as.matrix(newdata)
    pcac <- object@pcaobj[[1]]@center
    if(length(pcac)>0 & ncol(x) != length(pcac))
        stop("wrong number of variables")

    tt <- .testmodel(object, x)
    assgn <- .assigngroup(tt$odsc, tt$sdsc, method=method, gamma=gamma)
    cl <- factor(levels(object@grp)[assgn], levels=levels(object@grp))

    colnames(tt$odsc) <- colnames(tt$sdsc) <- names(object@prior)
    ret <- new("PredictSimca", classification=cl, odsc=tt$odsc, sdsc=tt$sdsc)
    if(ct){
        ctab <- rrcov::mtxconfusion(object@grp, cl)
        ret@ct <- ctab
    }

    ret
})

setMethod("show", "PredictSimca", function(object){

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

    invisible(object)
})

setMethod("summary", "Simca", function(object, ...){
    new("SummarySimca", simcaobj=object)
})

setMethod("show", "SummarySimca", function(object){

    cat("\nCall:\n")
    print(object@simcaobj@call)

    digits = max(3, getOption("digits") - 3)
    cat("\nPrior Probabilities of Groups:\n")
    print(object@simcaobj@prior)

    cat("\nPca objects for Groups:\n")
    for(i in 1:length(object@simcaobj@pcaobj))
    {
        ss <- summary(object@simcaobj@pcaobj[[i]])
        show(ss)
    }

    invisible(object)
})
